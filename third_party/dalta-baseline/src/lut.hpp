#ifndef LUT_HPP
#define LUT_HPP


#include <boost/timer/progress_display.hpp>

#include "bitOper.hpp"
#include "gurobi_c++.h"
#include "headers.h"


class LUT_t
{
private:
    int nIn;                                         // #LUT inputs
    int nOut;                                        // #LUT outputs
    int nAddr;                                       // #addresses
    int maxOutVal;                                   // maximum output value
    int nFree;                                       // #free variables
    int nBound;                                      // #bound variables
    int nOverl;                                      // #overlap variables
    int ymin;                                        // value domain [ymin, ymax]
    int ymax;
    int admmReInitTime;
    int admmStuckTol;
    unsigned seed;                                   // random seed
    double MIPTime;                                  // time limitation for MIP solver
    double MIPGap;                                   // |BP-BF|/|BP|, BP, best objective by relaxing integer->continous; BF, best objective by current integer solver
    std::vector <uint64_t> table;                    // accurate look-up table
    std::vector <double> prob;                       // prob[i], probability that the i-th pattern occurs
    std::vector < std::vector <bool> > truthTables;  // truthTables[i], accurate truth table for the i-th output
    std::vector <uint64_t> appTable;                 // approximate look-up table
    std::vector < std::vector <int> > bestFreeSets;  // bestFreeSets[i], best free set for the i-th output
    std::vector < std::vector <int> > bestBoundSets; // bestBoundSets[i], best bound set for the i-th output
    std::vector <int> bestFreeSetIds;                // bestFreeIds[i], best free set id for the i-th output
    std::vector <uint64_t> freeTable;                // free set table
    std::vector <uint64_t> boundTable;               // bound set table

    LUT_t & operator = (const LUT_t &);
    LUT_t(LUT_t &);
    LUT_t(LUT_t &&);

public:
    explicit LUT_t(std::string inpFile, int nFree, int admmReInitTime,
                   int admmStuckTol, unsigned seed, double MIPTime,
                   double MIPGap) {
        FILE * file = fopen(inpFile.c_str(), "r");
#ifdef DEBUG
        assert(file != nullptr);
#endif
        int addr, dat;
        table.reserve(1024);
        ymax = -INT_MAX;
        ymin = INT_MAX;
        while (fscanf(file, "%d %d", &addr, &dat) != -1) {
            table.emplace_back(dat);
#ifdef DEBUG
            assert(addr >= 0 && addr < (1 << 20));
            assert(dat >= 0 && dat < (1 << 20));
#endif
            ymax = std::max(ymax, dat);
            ymin = std::min(ymin, dat);
        }
        fclose(file);
        nIn = floor(log2(table.size()));
        nOut = floor(log2(ymax)) + 1;
        // nOut = floor(log2(ymax)) + 2;
        // nOut = 16;
        // std::cout << ymax << "," << nOut << std::endl;
        nAddr = 1 << nIn;
        maxOutVal = (1 << nOut) - 1;
        this->nFree = nFree;
        this->nBound = nIn - nFree;
        this->nOverl = nFree + nBound - nIn;
        this->admmReInitTime = admmReInitTime;
        this->admmStuckTol = admmStuckTol;
        this->seed = seed;
        this->MIPTime = MIPTime;
        this->MIPGap = MIPGap;
#ifdef DEBUG
        // std::cout << "accurate table = " << Vect2String(table) << std::endl;
        assert((1 << nIn) == static_cast <int> (table.size()));
        assert(nIn < 20);
        assert((1 << nOut) - 1 >= ymax);
        assert(nFree >= 1 && nFree < nIn);
        assert(nBound >= 1 && nBound < nIn);
        assert(nOverl >= 0 && nOverl < nIn);
        std::cout << "nIn = " << nIn << std::endl;
        std::cout << "nOut = " << nOut << std::endl;
        std::cout << "nFree = " << nFree << std::endl;
        std::cout << "nBound = " << nBound << std::endl;
        std::cout << "nOverl = " << nOverl << std::endl;
        std::cout << "value domain = [" << ymin << "," << ymax << "]" << std::endl;
#endif
        SetPattProb();
        SplitTruthTable();
    }

    ~LUT_t() {
    }

    inline unsigned GetSeed() {
        return seed;
    }

    void SetPattProb()
    {
        prob.resize(table.size());
        double val = 1 / static_cast <double> (prob.size());
        for (auto & slice: prob)
            slice = val;
#ifdef DEBUG
        // std::cout << "probability table = " << Vect2String(prob) << std::endl;
#endif
    }

    void SplitTruthTable()
    {
        truthTables.resize(nOut);
        for (int i = 0; i < nOut; ++i) {
            truthTables[i].resize(nAddr);
            for (int j = 0; j < nAddr; ++j)
                truthTables[i][j] = GetBit(table[j], i);
        }
#ifdef DEBUG
        // for (int i = nOut - 1; i >= 0; --i)
        //     std::cout << "truth table[" << i << "] = " << Vect2String(truthTables[i]) << std::endl;
#endif
    }

    void ApprDec(int startBit, int nTry, int nComb, const std::string & solver) {
        if (solver == "separate-ilp")
            SepApprDec(nComb, true);
        else if (solver == "separate-legacy")
            SepApprDec(nComb, false);
        else if (solver == "joint-ilp")
            IncApprDec(startBit, nTry, nComb, true);
        else if (solver == "joint-heuristic")
            IncApprDec(startBit, nTry, nComb, false);
        else
            throw std::invalid_argument("unknown solver: " + solver);
    }

    void GetComb(std::vector < std::vector <int> > & freeSets, std::vector <int> & freeSet, int startPos, int depth) {
        if (depth == nFree) {
            freeSets.emplace_back(freeSet);
            return;
        }
        for (int i = startPos; i < nIn; i++){
            freeSet[depth] = i;
            GetComb(freeSets, freeSet, i + 1, depth + 1);
        }
    }

    void EnumVarOrd(std::vector < std::vector <int> > & freeSets, std::vector < std::vector <int> > & boundSets) {
#ifdef DEBUG
        assert(nBound <= 10 && nFree <= 10);
#endif
        std::vector <int> freeSetTmp(nFree);
        freeSets.reserve(8192);
        GetComb(freeSets, freeSetTmp, 0, 0);
        GetBoundFromFree(freeSets, boundSets);
#ifdef DEBUG
        assert(freeSets.size() == boundSets.size());
        // for (int i = 0; i < static_cast <int> (freeSets.size()); ++i)
        //     std::cout << "freeSets[" << i << "] = " << Vect2String(freeSets[i]) << ", " << "boundSets[" << i << "] = " << Vect2String(boundSets[i]) << std::endl;
#endif
    }

    void DrawVarOrd(const int nComb, std::vector < std::vector <int> > & freeSets, std::vector < std::vector <int> > & boundSets) {
#ifdef DEBUG
        assert(nComb < 10000);
        assert(nBound >= 1 && nFree >= 1);
#endif
        freeSets.resize(nComb);
        boost::random::mt19937 randEng(seed);
        std::vector <int> vars(nIn);
        for (int j = 0; j < nIn; ++j)
            vars[j] = j;
        for (int i = 0; i < nComb; ++i) {
            // std::cout << Vect2String(vars) << std::endl;
            auto & freeSet = freeSets[i];
            freeSet.resize(nFree);
            for (int j = 0; j < nFree; ++j)
                freeSet[j] = vars[j];
            std::sort(freeSet.begin(), freeSet.end());
            std::shuffle(vars.begin(), vars.end(), randEng);
        }
        GetBoundFromFree(freeSets, boundSets);
#ifdef DEBUG
        assert(freeSets.size() == boundSets.size());
        // for (int i = 0; i < static_cast <int> (freeSets.size()); ++i)
        //     std::cout << "freeSets[" << i << "] = " << Vect2String(freeSets[i]) << ", " << "boundSets[" << i << "] = " << Vect2String(boundSets[i]) << std::endl;
#endif
    }

    // void NonDisjEnumVarOrd(std::vector < std::vector <int> > & freeSets, std::vector < std::vector <int> > & boundSets) {
    //     assert(nBound <= 10 && nFree <= 10);
    //     assert(nOverl >= 0 && nOverl < nIn);
    //     assert(nOverl == 1);
    //     std::vector <int> tempFreeSet(nFree);
    //     std::vector < std::vector <int> > tempFreeSets;
    //     tempFreeSets.reserve(8192);
    //     GetComb(tempFreeSets, tempFreeSet, 0, 0);
    //     freeSets.resize(tempFreeSets.size() * nFree);
    //     boundSets.resize(tempFreeSets.size() * nFree);
    //     int count = 0;
    //     for (int i = 0; i < static_cast <int> (tempFreeSets.size()); ++i) {
    //         for (int overlId = 0; overlId < nFree; ++overlId) {
    //             freeSets[count].assign(tempFreeSets[i].begin(), tempFreeSets[i].end());
    //             GetBoundFromFree(freeSets[count], boundSets[count], overlId);
    //             // std::cout << count << ", " << overlId << ", " << Vect2String(freeSets[count]) << ", " << Vect2String(boundSets[count]) << std::endl;
    //             ++count;
    //         }
    //     }
    //     assert(freeSets.size() == boundSets.size());
    //     // for (int i = 0; i < static_cast <int> (freeSets.size()); ++i)
    //     //     std::cout << "freeSets[" << i << "] = " << Vect2String(freeSets[i]) << ", " << "boundSets[" << i << "] = " << Vect2String(boundSets[i]) << std::endl;
    // }

    double SolveSepOpt(std::vector <bool> & truthTable, std::vector <int> & freeSet, std::vector <int> & boundSet, std::vector <bool> & bestRowPat, std::vector <int> & bestRowTyp) {
        // generate 2D table
        std::vector < std::vector <bool> > truTab2D;
        std::vector < std::vector <double> > prob2D;
        Gen2DTable(truthTable, freeSet, boundSet, truTab2D, prob2D);

        // linear programming
#ifdef DEBUG
        assert(nFree <= 10);
        assert(nBound <= 10);
#endif
        int m = 1 << nFree; // m rows
        int n = 1 << nBound; // n columns
        int nType = 3; // 0 ~ const, 1 ~ p, 2 ~ not p
        double ret = 0.0;
        GRBVar * pat = nullptr;
        GRBVar * diffRowErr = nullptr;
        GRBVar * rowErr = nullptr;
        GRBVar * pointDiffP = nullptr;
        GRBVar * pointErrP = nullptr;
        GRBVar * pointDiffPComp = nullptr;
        GRBVar * pointErrPComp = nullptr;
        try {
            // Create an environment
            GRBEnv env = GRBEnv(true);
            env.set(GRB_IntParam_LogToConsole, 0);
            env.start();

            // Create an empty model
            GRBModel model = GRBModel(env);
            model.set(GRB_DoubleParam_MIPGap, MIPGap);

            // Create variables
            pat = model.addVars(n);
            diffRowErr = model.addVars(m * nType);
            rowErr = model.addVars(m);
            pointDiffP = model.addVars(m * n);
            pointErrP = model.addVars(m * n);
            pointDiffPComp = model.addVars(m * n);
            pointErrPComp = model.addVars(m * n);
            for (int i = 0; i < n; ++i) {
                pat[i].set(GRB_DoubleAttr_LB, 0); // pattern vector
                pat[i].set(GRB_DoubleAttr_UB, 1);
                pat[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
            }
            for (int i = 0; i < m; ++i) {
                rowErr[i].set(GRB_DoubleAttr_LB, 0); // rowErr[i], minimum error of i-th row
                rowErr[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                rowErr[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                for (int j = 0; j < nType; ++j) {
                    diffRowErr[i * nType + j].set(GRB_DoubleAttr_LB, 0); // diffRowErr[i * nType + j], error of i-th row for j-th row type
                    diffRowErr[i * nType + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                    diffRowErr[i * nType + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                }
                for (int j = 0; j < n; ++j) {
                    pointDiffP[i * n + j].set(GRB_DoubleAttr_LB, -GRB_INFINITY); // pointDiffP[i * n + j], difference at point (i, j) for row type p
                    pointDiffP[i * n + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                    pointDiffP[i * n + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                    pointErrP[i * n + j].set(GRB_DoubleAttr_LB, 0);
                    pointErrP[i * n + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                    pointErrP[i * n + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                    pointDiffPComp[i * n + j].set(GRB_DoubleAttr_LB, -GRB_INFINITY); // pointDiffPComp[i * n + j], difference at point (i, j) for row type not p
                    pointDiffPComp[i * n + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                    pointDiffPComp[i * n + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                    pointErrPComp[i * n + j].set(GRB_DoubleAttr_LB, 0);
                    pointErrPComp[i * n + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                    pointErrPComp[i * n + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                }
            }

            // Set objective
            GRBLinExpr patErr = 0;
            for (int i = 0; i < m; ++i)
                patErr += rowErr[i];
            model.setObjective(patErr, GRB_MINIMIZE);

            // Add constraints
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    model.addConstr(pointDiffP[i * n + j] == truTab2D[i][j] - pat[j]);
                    model.addConstr(pointDiffPComp[i * n + j] == truTab2D[i][j] + pat[j] - 1);
                    model.addGenConstrAbs(pointErrP[i * n + j], pointDiffP[i * n + j]);
                    model.addGenConstrAbs(pointErrPComp[i * n + j], pointDiffPComp[i * n + j]);
                }
            }
            for (int i = 0; i < m; ++i) {
                // all 0/1
                double rowErr0 = 0;
                double rowErr1 = 0;
                for (int j = 0; j < n; ++j) {
                    rowErr0 += prob2D[i][j] * truTab2D[i][j];
                    rowErr1 += prob2D[i][j] * (1 - truTab2D[i][j]);
                }
                double minRowErr01 = std::min(rowErr0, rowErr1);
                model.addConstr(diffRowErr[i * nType] == minRowErr01);
                // p and p's complement
                GRBLinExpr rowErrP = 0;
                GRBLinExpr rowErrPComp = 0;
                for (int j = 0; j < n; ++j) {
                    rowErrP += prob2D[i][j] * pointErrP[i * n + j];
                    rowErrPComp += prob2D[i][j] * pointErrPComp[i * n + j];
                }
                model.addConstr(diffRowErr[i * nType + 1] == rowErrP);
                model.addConstr(diffRowErr[i * nType + 2] == rowErrPComp);
                model.addGenConstrMin(rowErr[i], diffRowErr + i * nType, nType);
            }

            // Optimize
            model.optimize();

            // Measure accurate error
            std::vector <int> rowPat(n, 0);
            for (int i = 0; i < n; ++i)
                rowPat[i] = static_cast <int> (round(pat[i].get(GRB_DoubleAttr_X)));
            std::vector <int> rowTyp(m, 0);
            double patErrVal = 0;
            for (int i = 0; i < m; ++i) {
                double rowErrVal = DBL_MAX;
                // try zero
                double tmpErr = 0;
                for (int j = 0; j < n; ++j)
                    tmpErr += prob2D[i][j] * truTab2D[i][j];
                if (tmpErr < rowErrVal) {
                    rowErrVal = tmpErr;
                    rowTyp[i] = 0;
                }
                // try one
                tmpErr = 0;
                for (int j = 0; j < n; ++j)
                    tmpErr += prob2D[i][j] * (1 - truTab2D[i][j]);
                if (tmpErr < rowErrVal) {
                    rowErrVal = tmpErr;
                    rowTyp[i] = 1;
                }
                // try pat
                tmpErr = 0;
                for (int j = 0; j < n; ++j)
                    tmpErr += prob2D[i][j] * abs(truTab2D[i][j] - rowPat[j]);
                if (tmpErr < rowErrVal) {
                    rowErrVal = tmpErr;
                    rowTyp[i] = 2;
                }
                // try 1 - pat
                tmpErr = 0;
                for (int j = 0; j < n; ++j)
                    tmpErr += prob2D[i][j] * abs(truTab2D[i][j] + rowPat[j] - 1);
                if (tmpErr < rowErrVal) {
                    rowErrVal = tmpErr;
                    rowTyp[i] = 3;
                }
                patErrVal += rowErrVal;
            }
            ret = patErrVal;

            // update
            bestRowPat.assign(rowPat.begin(), rowPat.end());
            bestRowTyp.assign(rowTyp.begin(), rowTyp.end());

        } catch(GRBException e) {
            std::cout << "Error code = " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }
        delete []pat;
        delete []rowErr;
        delete []diffRowErr;
        delete []pointDiffP;
        delete []pointErrP;
        delete []pointDiffPComp;
        delete []pointErrPComp;

        return ret;
    }

    // Algebraically equivalent formulation of DALTA's separate-mode ILP.
    // It avoids nested general MIN/ABS constraints, whose implementation in
    // recent Gurobi releases is dramatically slower than Gurobi 9.1 used by
    // the paper.  All four candidate row errors are linear because the exact
    // truth-table entries are constants and the pattern entries are binary.
    double SolveSepOptPortable(std::vector <bool> & truthTable,
                               std::vector <int> & freeSet,
                               std::vector <int> & boundSet,
                               std::vector <bool> & bestRowPat,
                               std::vector <int> & bestRowTyp) {
        std::vector < std::vector <bool> > truTab2D;
        std::vector < std::vector <double> > prob2D;
        Gen2DTable(truthTable, freeSet, boundSet, truTab2D, prob2D);

        const int nRow = 1 << nFree;
        const int nCol = 1 << nBound;
        double ret = 0.0;
        GRBVar * pat = nullptr;
        GRBVar * typ = nullptr;
        GRBVar * rowErr = nullptr;

        try {
            GRBEnv env = GRBEnv(true);
            env.set(GRB_IntParam_LogToConsole, 0);
            env.start();
            GRBModel model = GRBModel(env);
            model.set(GRB_DoubleParam_MIPGap, MIPGap);
            model.set(GRB_DoubleParam_TimeLimit, MIPTime);

            pat = model.addVars(nCol);
            typ = model.addVars(nRow * 4);
            rowErr = model.addVars(nRow);
            for (int j = 0; j < nCol; ++j)
                pat[j].set(GRB_CharAttr_VType, GRB_BINARY);
            for (int i = 0; i < nRow; ++i) {
                rowErr[i].set(GRB_DoubleAttr_LB, 0.0);
                rowErr[i].set(GRB_DoubleAttr_UB, 1.0);
                for (int d = 0; d < 4; ++d)
                    typ[i * 4 + d].set(GRB_CharAttr_VType, GRB_BINARY);
            }

            GRBLinExpr objective = 0.0;
            for (int i = 0; i < nRow; ++i) {
                objective += rowErr[i];
                model.addConstr(typ[i * 4] + typ[i * 4 + 1] +
                                typ[i * 4 + 2] + typ[i * 4 + 3] == 1.0);

                GRBLinExpr err0 = 0.0;
                GRBLinExpr err1 = 0.0;
                GRBLinExpr errPat = 0.0;
                GRBLinExpr errComp = 0.0;
                double bigM = 0.0;
                for (int j = 0; j < nCol; ++j) {
                    const double p = prob2D[i][j];
                    const bool exact = truTab2D[i][j];
                    bigM += p;
                    err0 += p * static_cast<double>(exact);
                    err1 += p * static_cast<double>(!exact);
                    errPat += p * (exact ? (1.0 - pat[j]) : pat[j]);
                    errComp += p * (exact ? pat[j] : (1.0 - pat[j]));
                }
                model.addConstr(rowErr[i] >= err0 - bigM * (1.0 - typ[i * 4]));
                model.addConstr(rowErr[i] >= err1 - bigM * (1.0 - typ[i * 4 + 1]));
                model.addConstr(rowErr[i] >= errPat - bigM * (1.0 - typ[i * 4 + 2]));
                model.addConstr(rowErr[i] >= errComp - bigM * (1.0 - typ[i * 4 + 3]));
            }
            // V and its complement describe the same decomposition after
            // swapping row types 2 and 3; remove that exact symmetry.
            model.addConstr(pat[0] == 0.0);
            model.setObjective(objective, GRB_MINIMIZE);
            model.optimize();

            bestRowPat.resize(nCol);
            for (int j = 0; j < nCol; ++j)
                bestRowPat[j] = pat[j].get(GRB_DoubleAttr_X) >= 0.5;
            bestRowTyp.resize(nRow);
            for (int i = 0; i < nRow; ++i) {
                bestRowTyp[i] = 0;
                for (int d = 1; d < 4; ++d)
                    if (typ[i * 4 + d].get(GRB_DoubleAttr_X) >
                        typ[i * 4 + bestRowTyp[i]].get(GRB_DoubleAttr_X))
                        bestRowTyp[i] = d;
                ret += rowErr[i].get(GRB_DoubleAttr_X);
            }
        } catch (GRBException & e) {
            std::cerr << "Gurobi error " << e.getErrorCode() << ": "
                      << e.getMessage() << std::endl;
            throw;
        }

        delete [] pat;
        delete [] typ;
        delete [] rowErr;
        return ret;
    }

    double SolveIncOpt(int tabId, std::vector <int> & freeSet, std::vector <int> & boundSet, std::vector <bool> & bestRowPat, std::vector <int> & bestRowTyp) {
#ifdef DEBUG
        assert(nFree <= 10);
        assert(nBound <= 10);
#endif
        // generate 2D table
        std::vector <bool> & truthTable = truthTables[tabId];
        std::vector < std::vector <bool> > truTab2D;
        std::vector < std::vector <double> > prob2D;
        std::vector < std::vector <int> > accTab2D;
        std::vector < std::vector <int> > appTab2D;
        Gen2DTable(truthTable, freeSet, boundSet, truTab2D, prob2D);
        Gen2DTable(table, freeSet, boundSet, accTab2D);
        Gen2DTable(appTable, freeSet, boundSet, appTab2D);
        for (auto & row: appTab2D) {
            for (auto & point: row)
                ResetBit(point, tabId);
        }

        // linear programming
        int m = 1 << nFree;
        int n = 1 << nBound;
        int nType = 3;
        int weight = 1 << tabId;
        double ret = 0.0;
        GRBVar * pat = nullptr;
        GRBVar * diffRowErr = nullptr;
        GRBVar * rowErr = nullptr;
        GRBVar * pointDiffP = nullptr;
        GRBVar * pointErrP = nullptr;
        GRBVar * pointDiffPComp = nullptr;
        GRBVar * pointErrPComp = nullptr;

        try {
            // create an environment
            GRBEnv env = GRBEnv(true);
            env.set(GRB_IntParam_LogToConsole, 0);
            env.start();

            // create an empty model
            GRBModel model = GRBModel(env);
            model.set(GRB_DoubleParam_MIPGap, MIPGap);

            // create variables
            pat = model.addVars(n);
            diffRowErr = model.addVars(m * nType);
            rowErr = model.addVars(m);
            pointDiffP = model.addVars(m * n);
            pointErrP = model.addVars(m * n);
            pointDiffPComp = model.addVars(m * n);
            pointErrPComp = model.addVars(m * n);
            for (int i = 0; i < n; ++i) {
                pat[i].set(GRB_DoubleAttr_LB, 0);
                pat[i].set(GRB_DoubleAttr_UB, 1);
                // pat[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                pat[i].set(GRB_CharAttr_VType, GRB_BINARY);
            }
            for (int i = 0; i < m; ++i) {
                rowErr[i].set(GRB_DoubleAttr_LB, 0);
                rowErr[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                rowErr[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                for (int j = 0; j < nType; ++j) {
                    diffRowErr[i * nType + j].set(GRB_DoubleAttr_LB, 0);
                    diffRowErr[i * nType + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                    diffRowErr[i * nType + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                }
                for (int j = 0; j < n; ++j) {
                    pointDiffP[i * n + j].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
                    pointDiffP[i * n + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                    pointDiffP[i * n + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                    pointErrP[i * n + j].set(GRB_DoubleAttr_LB, 0);
                    pointErrP[i * n + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                    pointErrP[i * n + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                    pointDiffPComp[i * n + j].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
                    pointDiffPComp[i * n + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                    pointDiffPComp[i * n + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                    pointErrPComp[i * n + j].set(GRB_DoubleAttr_LB, 0);
                    pointErrPComp[i * n + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                    pointErrPComp[i * n + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                }
            }

            // set objective
            GRBLinExpr patErr = 0;
            for (int i = 0; i < m; ++i)
                patErr += rowErr[i];
            model.setObjective(patErr, GRB_MINIMIZE);

            // add constraints
            // model.addConstr(pat[0] == 0); // from symmetry
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    model.addConstr(pointDiffP[i * n + j] == accTab2D[i][j] - (pat[j] * weight + appTab2D[i][j])); // ???
                    model.addConstr(pointDiffPComp[i * n + j] == accTab2D[i][j] -  ((1 - pat[j]) * weight + appTab2D[i][j])); // ???
                    model.addGenConstrAbs(pointErrP[i * n + j], pointDiffP[i * n + j]);
                    model.addGenConstrAbs(pointErrPComp[i * n + j], pointDiffPComp[i * n + j]);
                }
            }
            for (int i = 0; i < m; ++i) {
                // all 0/1
                double rowErr0 = 0;
                double rowErr1 = 0;
                for (int j = 0; j < n; ++j) {
                    rowErr0 += prob2D[i][j] * abs(accTab2D[i][j] - (appTab2D[i][j] + 0 * weight));
                    rowErr1 += prob2D[i][j] * abs(accTab2D[i][j] - (appTab2D[i][j] + 1 * weight));
                }
                double minRowErr01 = std::min(rowErr0, rowErr1);
                model.addConstr(diffRowErr[i * nType] == minRowErr01);
                // p and p's complement
                GRBLinExpr rowErrP = 0;
                GRBLinExpr rowErrPComp = 0;
                for (int j = 0; j < n; ++j) {
                    rowErrP += prob2D[i][j] * pointErrP[i * n + j];
                    rowErrPComp += prob2D[i][j] * pointErrPComp[i * n + j];
                }
                model.addConstr(diffRowErr[i * nType + 1] == rowErrP);
                model.addConstr(diffRowErr[i * nType + 2] == rowErrPComp);
                model.addGenConstrMin(rowErr[i], diffRowErr + i * nType, nType);
            }

            // optimize
            model.optimize();

            // measure accurate error
            std::vector <int> rowPat(n, 0);
            for (int i = 0; i < n; ++i) {
                rowPat[i] = static_cast <int> (round(pat[i].get(GRB_DoubleAttr_X)));
#ifdef DEBUG
                const double epsilon = 1e-5;
                if (abs(rowPat[i] - pat[i].get(GRB_DoubleAttr_X)) >= epsilon) {
                    printf("%d, %.9lf\n", rowPat[i], pat[i].get(GRB_DoubleAttr_X));
                    assert(0);
                }
#endif
            }

            // get row type vector
            std::vector <int> rowTyp(m, 0);
            double patErrVal = 0;
            for (int i = 0; i < m; ++i) {
                double rowErrVal = DBL_MAX;
                // try zero
                double tmpErr = 0;
                for (int j = 0; j < n; ++j)
                    tmpErr += prob2D[i][j] * abs(accTab2D[i][j] - (appTab2D[i][j] + 0 * weight));
                if (tmpErr < rowErrVal) {
                    rowErrVal = tmpErr;
                    rowTyp[i] = 0;
                }
                // try one
                tmpErr = 0;
                for (int j = 0; j < n; ++j)
                    tmpErr += prob2D[i][j] * abs(accTab2D[i][j] - (appTab2D[i][j] + 1 * weight));
                if (tmpErr < rowErrVal) {
                    rowErrVal = tmpErr;
                    rowTyp[i] = 1;
                }
                // try pat
                tmpErr = 0;
                for (int j = 0; j < n; ++j)
                    tmpErr += prob2D[i][j] * abs(accTab2D[i][j] -  (rowPat[j] * weight + appTab2D[i][j]));
                if (tmpErr < rowErrVal) {
                    rowErrVal = tmpErr;
                    rowTyp[i] = 2;
                }
                // try 1 - pat
                tmpErr = 0;
                for (int j = 0; j < n; ++j)
                    tmpErr += prob2D[i][j] * abs(accTab2D[i][j] -  ((1 - rowPat[j]) * weight + appTab2D[i][j]));
                if (tmpErr < rowErrVal) {
                    rowErrVal = tmpErr;
                    rowTyp[i] = 3;
                }
                patErrVal += rowErrVal;
                // std::cout << i << "," << rowErrVal << std::endl;
            }
            ret = patErrVal;

#ifdef DEBUG
            const double epsilon = 1e-5;
            // std::cout << patErrVal << "," << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
            assert(abs(patErrVal - model.get(GRB_DoubleAttr_ObjVal)) < epsilon);
#endif

            // update
            bestRowPat.assign(rowPat.begin(), rowPat.end());
            bestRowTyp.assign(rowTyp.begin(), rowTyp.end());

        } catch(GRBException e) {
            std::cout << "Error code = " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }
        delete []pat;
        delete []rowErr;
        delete []diffRowErr;
        delete []pointDiffP;
        delete []pointErrP;
        delete []pointDiffPComp;
        delete []pointErrPComp;

        // std::cout << "row pattern " << Vect2String(bestRowPat) << std::endl;
        // std::cout << "row type " << Vect2String(bestRowTyp) << std::endl;
        // std::cout << "objective = " << ret << std::endl;
        return ret;
    }

    // double NonDisjSolveIncOpt(const int tabId, const std::vector <int> & freeSet, const std::vector <int> & boundSet, std::vector <bool> & bestRowPat, std::vector <int> & bestRowTyp) {
    //     assert(tabId >= 0 && tabId < nOut);
    //     assert(static_cast <int> (freeSet.size()) == nFree);
    //     assert(static_cast <int> (boundSet.size()) == nBound);
    //     // generate 2D table
    // }

    int64_t GetTypeAndSEDFromPatt(const int tabId, const std::vector <std::vector <int> > & accTab2D, const std::vector <std::vector <int> > & appTab2D, const std::vector <bool> & pattVect, std::vector <int> & typeVect) {
#ifdef DEBUG
        assert(accTab2D.size());
        assert(accTab2D.size() == typeVect.size());
        assert(accTab2D[0].size());
        assert(accTab2D[0].size() == pattVect.size());
        assert(accTab2D.size() == appTab2D.size());
        assert(accTab2D[0].size() == appTab2D[0].size());
#endif
        int nRow = static_cast <int> (accTab2D.size());
        int nCol = static_cast <int> (accTab2D[0].size());
        int64_t sed = 0;
        boost::random::mt19937 randEng(seed);
        boost::uniform_int <int> unDistr(0, 1);
        for (int i = 0; i < nRow; ++i) {
            // try all 0
            int rowSed = 0;
            for (int j = 0; j < nCol; ++j)
                rowSed += abs(appTab2D[i][j] - accTab2D[i][j]);
            int bestRowSed = rowSed;
            typeVect[i] = 0;
            // try all 1
            rowSed = 0;
            for (int j = 0; j < nCol; ++j)
                rowSed += abs(appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]);
            if (rowSed < bestRowSed) {
                bestRowSed = rowSed;
                typeVect[i] = 1;
            }
            else if (rowSed == bestRowSed && unDistr(randEng))
                typeVect[i] = 1;
            // try pattVec
            rowSed = 0;
            for (int j = 0; j < nCol; ++j)
                rowSed += abs(appTab2D[i][j] + (pattVect[j] << tabId) - accTab2D[i][j]);
            if (rowSed < bestRowSed) {
                bestRowSed = rowSed;
                typeVect[i] = 2;
            }
            else if (rowSed == bestRowSed && unDistr(randEng))
                typeVect[i] = 2;
            // try 1 - pattVec
            rowSed = 0;
            for (int j = 0; j < nCol; ++j)
                rowSed += abs(appTab2D[i][j] + ((1 - pattVect[j]) << tabId) - accTab2D[i][j]);
            if (rowSed < bestRowSed) {
                bestRowSed = rowSed;
                typeVect[i] = 3;
            }
            else if (rowSed == bestRowSed && unDistr(randEng))
                typeVect[i] = 3;
            sed += bestRowSed;
        }
        return sed;
    }

//     int64_t GetTypeAndSSEFromPatt(const int tabId, const std::vector <std::vector <int> > & accTab2D, const std::vector <std::vector <int> > & appTab2D, const std::vector <bool> & pattVect, std::vector <int> & typeVect) {
// #ifdef DEBUG
//         assert(accTab2D.size());
//         assert(accTab2D.size() == typeVect.size());
//         assert(accTab2D[0].size());
//         assert(accTab2D[0].size() == pattVect.size());
//         assert(accTab2D.size() == appTab2D.size());
//         assert(accTab2D[0].size() == appTab2D[0].size());
// #endif
//         int nRow = static_cast <int> (accTab2D.size());
//         int nCol = static_cast <int> (accTab2D[0].size());
//         int64_t sed = 0;
//         boost::random::mt19937 randEng(seed);
//         boost::uniform_int <int> unDistr(0, 1);
//         for (int i = 0; i < nRow; ++i) {
//             // try all 0
//             int64_t rowSed = 0;
//             for (int j = 0; j < nCol; ++j)
//                 // rowSed += abs(appTab2D[i][j] - accTab2D[i][j]);
//                 rowSed += (appTab2D[i][j] - accTab2D[i][j]) * (appTab2D[i][j] - accTab2D[i][j]);
//             int64_t bestRowSed = rowSed;
//             typeVect[i] = 0;
//             // try all 1
//             rowSed = 0;
//             for (int j = 0; j < nCol; ++j)
//                 // rowSed += abs(appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]);
//                 rowSed += (appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]) * (appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]);
//             if (rowSed < bestRowSed) {
//                 bestRowSed = rowSed;
//                 typeVect[i] = 1;
//             }
//             else if (rowSed == bestRowSed && unDistr(randEng))
//                 typeVect[i] = 1;
//             // try pattVec
//             rowSed = 0;
//             for (int j = 0; j < nCol; ++j)
//                 // rowSed += abs(appTab2D[i][j] + (pattVect[j] << tabId) - accTab2D[i][j]);
//                 rowSed += (appTab2D[i][j] + (pattVect[j] << tabId) - accTab2D[i][j]) * (appTab2D[i][j] + (pattVect[j] << tabId) - accTab2D[i][j]);
//             if (rowSed < bestRowSed) {
//                 bestRowSed = rowSed;
//                 typeVect[i] = 2;
//             }
//             else if (rowSed == bestRowSed && unDistr(randEng))
//                 typeVect[i] = 2;
//             // try 1 - pattVec
//             rowSed = 0;
//             for (int j = 0; j < nCol; ++j)
//                 // rowSed += abs(appTab2D[i][j] + ((1 - pattVect[j]) << tabId) - accTab2D[i][j]);
//                 rowSed += (appTab2D[i][j] + ((1 - pattVect[j]) << tabId) - accTab2D[i][j]) * (appTab2D[i][j] + ((1 - pattVect[j]) << tabId) - accTab2D[i][j]);
//             if (rowSed < bestRowSed) {
//                 bestRowSed = rowSed;
//                 typeVect[i] = 3;
//             }
//             else if (rowSed == bestRowSed && unDistr(randEng))
//                 typeVect[i] = 3;
//             sed += bestRowSed;
//         }
//         return sed;
//     }

    int64_t GetPattAndSEDFromType(const int tabId, const std::vector <std::vector <int> > & accTab2D, const std::vector <std::vector <int> > & appTab2D, std::vector <bool> & pattVect, const std::vector <int> & typeVect) {
#ifdef DEBUG
        assert(accTab2D.size());
        assert(accTab2D.size() == typeVect.size());
        assert(accTab2D[0].size());
        assert(accTab2D[0].size() == pattVect.size());
        assert(accTab2D.size() == appTab2D.size());
        assert(accTab2D[0].size() == appTab2D[0].size());
#endif
        int nRow = static_cast <int> (accTab2D.size());
        int nCol = static_cast <int> (accTab2D[0].size());
        int64_t sed = 0;
        boost::random::mt19937 randEng(seed);
        boost::uniform_int <int> unDistr(0, 1);
        for (int j = 0; j < nCol; ++j) {
            // try 0
            int colSed = 0;
            for (int i = 0; i < nRow; ++i) {
                if (typeVect[i] == 0 || typeVect[i] == 2) // 0, patt
                    colSed += abs(appTab2D[i][j] - accTab2D[i][j]);
                else if (typeVect[i] == 1 || typeVect[i] == 3) // 1, 1 - patt
                    colSed += abs(appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]);
                else
                    assert(0);
            }
            int bestColSed = colSed;
            pattVect[j] = 0;
            // try 1
            colSed = 0;
            for (int i = 0; i < nRow; ++i) {
                if (typeVect[i] == 0 || typeVect[i] == 3) // 0, 1 - patt
                    colSed += abs(appTab2D[i][j] - accTab2D[i][j]);
                else if (typeVect[i] == 1 || typeVect[i] == 2) // 1, patt
                    colSed += abs(appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]);
                else
                    assert(0);
            }
            if (colSed < bestColSed) {
                bestColSed = colSed;
                pattVect[j] = 1;
            }
            else if (colSed == bestColSed && unDistr(randEng))
                pattVect[j] = 1;
            sed += bestColSed;
        }
        return sed;
    }

//     int64_t GetPattAndSSEFromType(const int tabId, const std::vector <std::vector <int> > & accTab2D, const std::vector <std::vector <int> > & appTab2D, std::vector <bool> & pattVect, const std::vector <int> & typeVect) {
// #ifdef DEBUG
//         assert(accTab2D.size());
//         assert(accTab2D.size() == typeVect.size());
//         assert(accTab2D[0].size());
//         assert(accTab2D[0].size() == pattVect.size());
//         assert(accTab2D.size() == appTab2D.size());
//         assert(accTab2D[0].size() == appTab2D[0].size());
// #endif
//         int nRow = static_cast <int> (accTab2D.size());
//         int nCol = static_cast <int> (accTab2D[0].size());
//         int64_t sed = 0;
//         boost::random::mt19937 randEng(seed);
//         boost::uniform_int <int> unDistr(0, 1);
//         for (int j = 0; j < nCol; ++j) {
//             // try 0
//             int64_t colSed = 0;
//             for (int i = 0; i < nRow; ++i) {
//                 if (typeVect[i] == 0 || typeVect[i] == 2) // 0, patt
//                     // colSed += abs(appTab2D[i][j] - accTab2D[i][j]);
//                     colSed += (appTab2D[i][j] - accTab2D[i][j]) * (appTab2D[i][j] - accTab2D[i][j]);
//                 else if (typeVect[i] == 1 || typeVect[i] == 3) // 1, 1 - patt
//                     // colSed += abs(appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]);
//                     colSed += (appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]) * (appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]);
//                 else
//                     assert(0);
//             }
//             int64_t bestColSed = colSed;
//             pattVect[j] = 0;
//             // try 1
//             colSed = 0;
//             for (int i = 0; i < nRow; ++i) {
//                 if (typeVect[i] == 0 || typeVect[i] == 3) // 0, 1 - patt
//                     // colSed += abs(appTab2D[i][j] - accTab2D[i][j]);
//                     colSed += (appTab2D[i][j] - accTab2D[i][j]) * (appTab2D[i][j] - accTab2D[i][j]);
//                 else if (typeVect[i] == 1 || typeVect[i] == 2) // 1, patt
//                     // colSed += abs(appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]);
//                     colSed += (appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]) * (appTab2D[i][j] + (1 << tabId) - accTab2D[i][j]);
//                 else
//                     assert(0);
//             }
//             if (colSed < bestColSed) {
//                 bestColSed = colSed;
//                 pattVect[j] = 1;
//             }
//             else if (colSed == bestColSed && unDistr(randEng))
//                 pattVect[j] = 1;
//             sed += bestColSed;
//         }
//         return sed;
//     }

    double GetMedFromPattAndTypeVect(const std::vector <std::vector <bool> > & truTab2D, std::vector <bool> & pattVect, const std::vector <int> & typeVect) {
#ifdef DEBUG
        assert(truTab2D.size());
        assert(truTab2D.size() == typeVect.size());
        assert(truTab2D[0].size());
        assert(truTab2D[0].size() == pattVect.size());
#endif
        int nRow = static_cast <int> (truTab2D.size());
        int nCol = static_cast <int> (truTab2D[0].size());
        double med = 0.0;
        for (int i = 0; i < nRow; ++i) {
            if (typeVect[i] == 0) {
                for (int j = 0; j < nCol; ++j)
                    med += truTab2D[i][j];
            }
            else if (typeVect[i] == 1) {
                for (int j = 0; j < nCol; ++j)
                    med += (1 - truTab2D[i][j]);
            }
            else if (typeVect[i] == 2) {
                for (int j = 0; j < nCol; ++j)
                    med += abs(pattVect[j] - truTab2D[i][j]);
            }
            else if (typeVect[i] == 3) {
                for (int j = 0; j < nCol; ++j)
                    med += abs(1 - pattVect[j] - truTab2D[i][j]);
            }
        }
        med /= (nRow * nCol);
        return med;
    }

    double SolveIncOptADMM(const int tabId, const std::vector <int> & freeSet, const std::vector <int> & boundSet, std::vector <bool> & bestRowPat, std::vector <int> & bestRowTyp) {
#ifdef DEBUG
        assert(nFree <= 10);
        assert(nBound <= 10);
#endif
        // generate 2D table
        std::vector < std::vector <int> > accTab2D;
        std::vector < std::vector <int> > appTab2D;
        Gen2DTable(table, freeSet, boundSet, accTab2D);
        Gen2DTable(appTable, freeSet, boundSet, appTab2D);
        for (auto & row: appTab2D) {
            for (auto & point: row)
                ResetBit(point, tabId);
        }

        // init
        int nRow = 1 << nFree;
        int nCol = 1 << nBound;
        std::vector <bool> pattVect(nCol, 0);
        std::vector <int> typeVect(nRow, 0);
        boost::random::mt19937 randEng(seed);
        boost::uniform_int <int> unDistr(0, 1);

        // admm
        int64_t bestSed = INT64_MAX;
        for (int round = 0; round < admmReInitTime; ++round) {
            for (int i = 0; i < nCol; ++i)
                pattVect[i] = unDistr(randEng);
            int64_t sed = GetTypeAndSEDFromPatt(tabId, accTab2D, appTab2D, pattVect, typeVect);
            // int64_t sed = GetTypeAndSSEFromPatt(tabId, accTab2D, appTab2D, pattVect, typeVect);
            int64_t oldSed = sed;
            // std::cout << "pattVect = " << Vect2String(pattVect, 0, ",") << ", typeVect = " << Vect2String(typeVect, 0, ",") << ", med = " << static_cast <double> (sed) / nRow / nCol << std::endl;
            int isOpt = 1; // 0, terminate; 1, optimize pattVec; 2, optimize typeVect
            int stuckCount = 0;
            while (isOpt) {
                if (isOpt == 1) {
                    sed = GetPattAndSEDFromType(tabId, accTab2D, appTab2D, pattVect, typeVect);
                    // sed = GetPattAndSSEFromType(tabId, accTab2D, appTab2D, pattVect, typeVect);
                    isOpt = 2;
                }
                else if (isOpt == 2) {
                    sed = GetTypeAndSEDFromPatt(tabId, accTab2D, appTab2D, pattVect, typeVect);
                    // sed = GetTypeAndSSEFromPatt(tabId, accTab2D, appTab2D, pattVect, typeVect);
                    isOpt = 1;
                }
                if (oldSed < sed)
                    assert(0);
                else if (oldSed == sed) {
                    stuckCount += 1;
                    if (stuckCount >= admmStuckTol)
                        isOpt = 0;
                }
                else {
                    oldSed = sed;
                    stuckCount = 0;
                }
            }
            // std::cout << "pattVect = " << Vect2String(pattVect, 0, ",") << ", typeVect = " << Vect2String(typeVect, 0, ",") << ", med = " << static_cast <double> (sed) / nRow / nCol << ", bettSolTime = " << bettSolTime << std::endl;
            if (sed < bestSed) {
                bestRowPat.assign(pattVect.begin(), pattVect.end());
                bestRowTyp.assign(typeVect.begin(), typeVect.end());
                bestSed = sed;
            }
            else if (sed == bestSed) {
                if (unDistr(randEng)) {
                    bestRowPat.assign(pattVect.begin(), pattVect.end());
                    bestRowTyp.assign(typeVect.begin(), typeVect.end());
                }
            }
        }
        // std::cout << "best pattVect = " << Vect2String(bestRowPat, 0, ",") << ", best typeVect = " << Vect2String(bestRowTyp, 0, ",") << ", med = " << static_cast <double> (bestSed) / nRow / nCol << std::endl;

        return static_cast <double> (bestSed) / nRow / nCol;
    }

    double SolveFineGrained(std::vector <int> & freeSet, std::vector <int> & boundSet, std::vector <bool> & bestRowPat, std::vector <int> & bestRowTyp) {
        assert(nBound <= 10 && nFree <= 10);
        assert(maxOutVal == (1 << nOut) - 1);

        // generate 2D table
        std::vector < std::vector <int> > valTab2D;
        std::vector < std::vector <double> > prob2D;
        Gen2DValTable(freeSet, boundSet, valTab2D, prob2D);

        // integer programming
        int m = 1 << nFree;
        int n = 1 << nBound;
        assert(static_cast <int> (bestRowPat.size()) == nOut * n);
        assert(static_cast <int> (bestRowTyp.size()) == nOut * m);
        double minErr = 0.0;
        GRBVar * pat = nullptr;
        GRBVar * typ = nullptr;
        GRBVar * ed = nullptr;
        GRBVar * diff = nullptr;

        try {
            // Create an environment
            GRBEnv env = GRBEnv(true);
            env.set(GRB_IntParam_LogToConsole, 1);
            env.set(GRB_DoubleParam_TimeLimit, MIPTime);
            env.start();

            // Create an empty model
            GRBModel model = GRBModel(env);
            model.set(GRB_DoubleParam_MIPGap, MIPGap);

            // Create variables
            pat = model.addVars(nOut * n);
            typ = model.addVars(nOut * m * 4);
            ed = model.addVars(m * n);
            diff = model.addVars(m * n);
            char varName[100];
            for (int i = 0; i < nOut * n; ++i) {
                pat[i].set(GRB_DoubleAttr_LB, 0);
                pat[i].set(GRB_DoubleAttr_UB, 1);
                pat[i].set(GRB_CharAttr_VType, GRB_BINARY);
                sprintf(varName, "c[%d][%d]", i / n, i % n);
                pat[i].set(GRB_StringAttr_VarName, std::string(varName));
            }
            for (int i = 0; i < nOut * m * 4; ++i) {
                typ[i].set(GRB_DoubleAttr_LB, 0);
                typ[i].set(GRB_DoubleAttr_UB, 1);
                typ[i].set(GRB_CharAttr_VType, GRB_BINARY);
                sprintf(varName, "s[%d][%d][%d]", i / (m * 4), i % (m * 4) / 4, i % 4);
                typ[i].set(GRB_StringAttr_VarName, std::string(varName));
            }
            for (int i = 0; i < m * n; ++i) {
                ed[i].set(GRB_DoubleAttr_LB, 0);
                ed[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                ed[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                sprintf(varName, "ed[%d][%d]", i / n, i % n);
                ed[i].set(GRB_StringAttr_VarName, std::string(varName));
                diff[i].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
                diff[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                diff[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                sprintf(varName, "diff[%d][%d]", i / n, i % n);
                diff[i].set(GRB_StringAttr_VarName, std::string(varName));
            }

            // Set objective
            GRBLinExpr nmed = 0;
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                    nmed += prob2D[i][j] * ed[i * n + j];
            model.setObjective(nmed, GRB_MINIMIZE);

            // Add constraints
            for (int k = 0; k < nOut; ++k) {
                for (int i = 0; i < m; ++i) {
                    int base = k * m * 4 + i * 4;
                    model.addConstr(typ[base] + typ[base + 1] + typ[base + 2] + typ[base + 3] == 1);
                }
            }
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    double accOut = valTab2D[i][j];
                    GRBQuadExpr appOut = 0;
                    for (int k = 0; k < nOut; ++k) {
                        int base = k * m * 4 + i * 4;
                        appOut += (1 << k) * (typ[base + 1] + typ[base + 2] * pat[k * n + j] + typ[base + 3] * (1 - pat[k * n + j]));
                    }
                    model.addQConstr(diff[i * n + j] == accOut - appOut);
                    model.addGenConstrAbs(ed[i * n + j], diff[i * n + j]);
                }
            }

            // Optimize
            model.optimize();

            for (int i = 0; i < nOut * n; ++i)
                bestRowPat[i] = static_cast <int> (round(pat[i].get(GRB_DoubleAttr_X)));
            for (int i = 0; i < nOut * m; ++i) {
                if (typ[i * 4].get(GRB_DoubleAttr_X))
                    bestRowTyp[i] = 0;
                else if (typ[i * 4 + 1].get(GRB_DoubleAttr_X))
                    bestRowTyp[i] = 1;
                else if (typ[i * 4 + 2].get(GRB_DoubleAttr_X))
                    bestRowTyp[i] = 2;
                else if (typ[i * 4 + 3].get(GRB_DoubleAttr_X))
                    bestRowTyp[i] = 3;
                else
                    assert(0);
            }
            minErr = model.get(GRB_DoubleAttr_ObjVal) / maxOutVal;
            std::cout << "best row pattern for the free/bound set: " << Vect2String(bestRowPat, 0, ",") << std::endl;
            std::cout << "best row type for the free/bound set: " << Vect2String(bestRowTyp, 0, ",") << std::endl;
            std::cout << "min error for the free/bound set: " << minErr << std::endl;
        } catch(GRBException e) {
            std::cout << "Error code = " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }

        delete []pat;
        delete []typ;
        delete []ed;
        delete []diff;

        return minErr;
    }

    double SolveFineGrainedPro(int startBit, std::vector < std::vector <int> > & freeSets, std::vector < std::vector <int> > & boundSets, std::vector <bool> & bestRowPat, std::vector <int> & bestRowTyp)  {
#ifdef DEBUG
        assert(startBit + 1 == static_cast <int> (freeSets.size()));
        assert(freeSets.size() == boundSets.size());
        assert(startBit >= 0 && startBit <= nOut - 1);
        assert(nBound <= 10 && nFree <= 10);
        assert(maxOutVal + 1 == 1 << nOut);
#endif

        int nRow = 1 << nFree; // # rows
        int nCol = 1 << nBound; // # columns
        assert(static_cast <int> (bestRowPat.size()) == (startBit + 1) * nCol);
        assert(static_cast <int> (bestRowTyp.size()) == (startBit + 1) * nRow);
        double minErr = 0.0;
        GRBVar * pat = nullptr;
        GRBVar * typ = nullptr;
        GRBVar * ed = nullptr;
        GRBVar * diff = nullptr;

        try {
            // create an environment
            GRBEnv env = GRBEnv(true);
            env.set(GRB_IntParam_LogToConsole, 1);
            env.set(GRB_DoubleParam_TimeLimit, MIPTime);
            env.start();

            // create an empty model
            GRBModel model = GRBModel(env);
            model.set(GRB_DoubleParam_MIPGap, MIPGap);

            // create variables
            pat = model.addVars((startBit + 1) * nCol);
            typ = model.addVars((startBit + 1) * nRow * 4);
            ed = model.addVars(nRow * nCol);
            diff = model.addVars(nRow * nCol);
            char varName[100];
            for (int i = 0; i < (startBit + 1) * nCol; ++i) {
                pat[i].set(GRB_DoubleAttr_LB, 0);
                pat[i].set(GRB_DoubleAttr_UB, 1);
                pat[i].set(GRB_CharAttr_VType, GRB_BINARY);
                sprintf(varName, "c[%d][%d]", i / nCol, i % nCol);
                pat[i].set(GRB_StringAttr_VarName, std::string(varName));
            }
            for (int i = 0; i < (startBit + 1) * nRow * 4; ++i) {
                typ[i].set(GRB_DoubleAttr_LB, 0);
                typ[i].set(GRB_DoubleAttr_UB, 1);
                typ[i].set(GRB_CharAttr_VType, GRB_BINARY);
                sprintf(varName, "s[%d][%d][%d]", i / (nRow * 4), i % (nRow * 4) / 4, i % 4);
                typ[i].set(GRB_StringAttr_VarName, std::string(varName));
            }
            assert(nRow * nCol == nAddr);
            for (int i = 0; i < nAddr; ++i) {
                ed[i].set(GRB_DoubleAttr_LB, 0);
                ed[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                ed[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                sprintf(varName, "ed[%d]", i);
                ed[i].set(GRB_StringAttr_VarName, std::string(varName));
                diff[i].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
                diff[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                diff[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                sprintf(varName, "diff[%d]", i);
                diff[i].set(GRB_StringAttr_VarName, std::string(varName));
            }

            // set objective
            GRBLinExpr med = 0;
            for (int i = 0; i < nAddr; ++i)
                med += ed[i];
            med = med / nAddr;
            model.setObjective(med, GRB_MINIMIZE);

            // add constraints
            for (int tabId  = 0; tabId <= startBit; ++tabId) {
                for (int row = 0; row < nRow; ++row) {
                    int base = tabId * nRow * 4 + row * 4;
                    model.addConstr(typ[base] + typ[base + 1] + typ[base + 2] + typ[base + 3] == 1);
                }
            }
            int mask = (1 << (startBit + 1)) - 1;
            std::vector <bool> addrBits(nIn, 0);
            for (int addr = 0; addr < nAddr; ++addr) {
                int accOut = table[addr] & mask;
                for (int bit = 0; bit < nIn; ++bit)
                    addrBits[bit] = GetBit(addr, bit);
                // std::cout << addr << "," << Vect2String(addrBits, 1, "") << std::endl;
                GRBQuadExpr appOut = 0;
                for (int tabId = 0; tabId <= startBit; ++tabId) {
                    auto & freeSet = freeSets[tabId];
                    auto & boundSet = boundSets[tabId];
                    assert(static_cast <int> (freeSet.size() + boundSet.size()) == nIn);
                    int row = 0;
                    int col = 0;
                    for (int varId = freeSet.size() - 1; varId >= 0; --varId) {
                        row <<= 1;
                        row += addrBits[freeSet[varId]];
                    }
                    for (int varId = boundSet.size() - 1; varId >= 0; --varId) {
                        col <<= 1;
                        col += addrBits[boundSet[varId]];
                    }
                    // std::cout << tabId << "," << Vect2String(freeSet, 0, ",") << "," << Vect2String(boundSet, 0, ",") << "," << row << "," << col << std::endl;
                    int base = tabId * nRow * 4 + row * 4;
                    appOut += (1 << tabId) * (typ[base + 1] + typ[base + 2] * pat[tabId * nCol + col] + typ[base + 3] * (1 - pat[tabId * nCol + col]));
                }
                model.addQConstr(diff[addr] == accOut - appOut);
                model.addGenConstrAbs(ed[addr], diff[addr]);
            }

            // write model
            model.write("test.lp");

            // optimize
            model.optimize();

            // output
            for (int i = 0; i < (startBit + 1) * nCol; ++i)
                bestRowPat[i] = static_cast <int> (round(pat[i].get(GRB_DoubleAttr_X)));
            for (int i = 0; i < (startBit + 1) * nRow; ++i) {
                if (typ[i * 4].get(GRB_DoubleAttr_X))
                    bestRowTyp[i] = 0;
                else if (typ[i * 4 + 1].get(GRB_DoubleAttr_X))
                    bestRowTyp[i] = 1;
                else if (typ[i * 4 + 2].get(GRB_DoubleAttr_X))
                    bestRowTyp[i] = 2;
                else if (typ[i * 4 + 3].get(GRB_DoubleAttr_X))
                    bestRowTyp[i] = 3;
                else
                    assert(0);
            }
            std::cout << "row pattern " << Vect2String(bestRowPat) << std::endl;
            std::cout << "row type " << Vect2String(bestRowTyp) << std::endl;
            std::cout << "objective = " << model.get(GRB_DoubleAttr_ObjVal) / maxOutVal << std::endl;
        } catch(GRBException e) {
            std::cout << "Error code = " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }

        delete []pat;
        delete []typ;
        delete []ed;
        delete []diff;

        return minErr;
    }

    void Gen2DTable(std::vector <bool> & truthTable, std::vector <int> & freeSet, std::vector <int> & boundSet, std::vector < std::vector <bool> > & truTab2D, std::vector < std::vector <double> > & prob2D) {
        int nRows = 1 << freeSet.size();
        int nCols = 1 << boundSet.size();
        truTab2D.resize(nRows);
        prob2D.resize(nRows);
        for (int i = 0; i < nRows; ++i) {
            truTab2D[i].resize(nCols);
            prob2D[i].resize(nCols);
            int base = 0;
            for (int k = 0; k < nFree; ++k)
                base += (GetBit(static_cast <uint32_t> (i), static_cast <uint32_t> (k)) << freeSet[k]);
            for (int j = 0; j < nCols; ++j) {
                int pos = base;
                for (int k = 0; k < nBound; ++k)
                    pos += (GetBit(static_cast <uint32_t> (j), static_cast <uint32_t> (k)) << boundSet[k]);
                truTab2D[i][j] = truthTable[pos];
                prob2D[i][j] = prob[pos];
            }
        }
    }

    void Gen2DTable(std::vector <bool> & truthTable, std::vector <int> & freeSet, std::vector <int> & boundSet, std::vector < std::vector <bool> > & truTab2D) {
        int nRows = 1 << freeSet.size();
        int nCols = 1 << boundSet.size();
        truTab2D.resize(nRows);
        for (int i = 0; i < nRows; ++i) {
            truTab2D[i].resize(nCols);
            int base = 0;
            for (int k = 0; k < nFree; ++k)
                base += (GetBit(static_cast <uint32_t> (i), static_cast <uint32_t> (k)) << freeSet[k]);
            for (int j = 0; j < nCols; ++j) {
                int pos = base;
                for (int k = 0; k < nBound; ++k)
                    pos += (GetBit(static_cast <uint32_t> (j), static_cast <uint32_t> (k)) << boundSet[k]);
                truTab2D[i][j] = truthTable[pos];
            }
        }
    }

    void Gen2DTable(const std::vector <uint64_t> & valTable, const std::vector <int> & freeSet, const std::vector <int> & boundSet, std::vector < std::vector <int> > & valTab2D) {
        int nRows = 1 << freeSet.size();
        int nCols = 1 << boundSet.size();
        valTab2D.resize(nRows);
        for (int i = 0; i < nRows; ++i) {
            valTab2D[i].resize(nCols);
            int base = 0;
            for (int k = 0; k < nFree; ++k)
                base += (GetBit(static_cast <uint32_t> (i), static_cast <uint32_t> (k)) << freeSet[k]);
            for (int j = 0; j < nCols; ++j) {
                int pos = base;
                for (int k = 0; k < nBound; ++k)
                    pos += (GetBit(static_cast <uint32_t> (j), static_cast <uint32_t> (k)) << boundSet[k]);
                valTab2D[i][j] = valTable[pos];
            }
        }
    }

    void Gen2DValTable(std::vector <int> & freeSet, std::vector <int> & boundSet, std::vector < std::vector <int> > & valTab2D, std::vector < std::vector <double> > & prob2D) {
        int nRows = 1 << freeSet.size();
        int nCols = 1 << boundSet.size();
        valTab2D.resize(nRows);
        prob2D.resize(nRows);
        for (int i = 0; i < nRows; ++i) {
            valTab2D[i].resize(nCols);
            prob2D[i].resize(nCols);
            int base = 0;
            for (int k = 0; k < nFree; ++k)
                base += (GetBit(i, k) << freeSet[k]);
            for (int j = 0; j < nCols; ++j) {
                int pos = base;
                for (int k = 0; k < nBound; ++k)
                    pos += (GetBit(j, k) << boundSet[k]);
                valTab2D[i][j] = table[pos];
                prob2D[i][j] = prob[pos];
            }
        }
    }

    void JoinAppTruTables(std::vector <bool> & bestRowPat, std::vector <int> & bestRowTyp, int tabId) {
#ifdef DEBUG
        assert(static_cast <int> (bestRowTyp.size()) == 1 << nFree);
        assert(static_cast <int> (bestRowPat.size()) == 1 << nBound);
#endif
        for (int i = 0; i < static_cast <int> (bestRowTyp.size()); ++i) {
            int base = 0;
            for (int k = 0; k < nFree; ++k)
                base += (GetBit(static_cast <uint32_t> (i), static_cast <uint32_t> (k)) << bestFreeSets[tabId][k]);
            for (int j = 0; j < static_cast <int> (bestRowPat.size()); ++j) {
                bool val = 0;
                switch (bestRowTyp[i]) {
                    case 0:
                        val = 0;
                        break;
                    case 1:
                        val = 1;
                        break;
                    case 2:
                        val = bestRowPat[j];
                        break;
                    case 3:
                        val = !bestRowPat[j];
                        break;
                    default:
                        assert(0);
                }
                int pos = base;
                for (int k = 0; k < nBound; ++k)
                    pos += (GetBit(static_cast <uint32_t> (j), static_cast <uint32_t> (k)) << bestBoundSets[tabId][k]);
                // if (i == 0 && j == 6) {
                //     std::cout << val << std::endl;
                //     std::cout << pos << std::endl;
                //     std::cout << appTable[pos] << std::endl;
                // }
                if (val)
                    SetBit(appTable[pos], tabId);
                else
                    ResetBit(appTable[pos], tabId);
                // if (i == 0 && j == 6)
                //     std::cout << appTable[pos] << std::endl;
            }
        }
    }

    void GenFreeBoundTable(std::vector <bool> & bestRowPat, std::vector <int> & bestRowTyp, int tabId) {
#ifdef DEBUG
        assert(boundTable.size() == bestRowPat.size());
#endif
        for (int i = 0; i < static_cast <int> (bestRowPat.size()); ++i) {
            if (bestRowPat[i])
                SetBit(boundTable[i], tabId);
            else
                ResetBit(boundTable[i], tabId);
        }
#ifdef DEBUG
        assert(freeTable.size() == (bestRowTyp.size() << 1));
#endif
        for (int i = 0; i < static_cast <int> (bestRowTyp.size()); ++i) {
            switch (bestRowTyp[i]) {
                case 0:
                    ResetBit(freeTable[i << 1], tabId);
                    ResetBit(freeTable[(i << 1) + 1], tabId);
                    break;
                case 1:
                    SetBit(freeTable[i << 1], tabId);
                    SetBit(freeTable[(i << 1) + 1], tabId);
                    break;
                case 2:
                    ResetBit(freeTable[i << 1], tabId);
                    SetBit(freeTable[(i << 1) + 1], tabId);
                    break;
                case 3:
                    SetBit(freeTable[i << 1], tabId);
                    ResetBit(freeTable[(i << 1) + 1], tabId);
                    break;
                default:
                    assert(0);
            }
        }
    }

    void SepApprDec(int nComb, bool portable) {
#ifdef DEBUG
        assert(nBound + nFree == nIn);
#endif
        // enumerate free sets and bound sets
        std::vector < std::vector <int> > freeSets;
        std::vector < std::vector <int> > boundSets;
        EnumVarOrd(freeSets, boundSets);
        if (nComb > 0 && nComb < static_cast<int>(freeSets.size())) {
            freeSets.resize(nComb);
            boundSets.resize(nComb);
        }

#ifdef DEBUG
        assert(nOut == static_cast <int> (truthTables.size()));
        assert(freeSets.size() < (1 << 20));
        assert(freeSets.size() == boundSets.size());
#endif
        // decompose each truth table
        appTable.resize(table.size(), 0);
        bestFreeSets.resize(nOut);
        bestBoundSets.resize(nOut);
        bestFreeSetIds.resize(nOut);
        freeTable.resize(1 << (nFree + 1), 0);
        boundTable.resize(1 << nBound, 0);
        for (int tabId = 0; tabId < nOut; ++tabId) {
#ifdef DEBUG
            std::cout << "**************** approximating " << tabId << "-th table ****************" << std::endl;
#endif
            double bestErr = DBL_MAX;
            std::vector <bool> bestRowPat;
            std::vector <int> bestRowTyp;
            std::vector <int> & bestFreeSet = bestFreeSets[tabId];
            std::vector <int> & bestBoundSet = bestBoundSets[tabId];
            // for each free/bound set pair, find the best decomposition
            for (int i = 0; i < static_cast <int> (freeSets.size()); ++i) {
                std::vector <bool> bestRowPatFixOrd(1 << nBound);
                std::vector <int> bestRowTypFixOrd((1 << nFree));
                double bestErrFixOrd = portable
                    ? SolveSepOptPortable(truthTables[tabId], freeSets[i], boundSets[i], bestRowPatFixOrd, bestRowTypFixOrd)
                    : SolveSepOpt(truthTables[tabId], freeSets[i], boundSets[i], bestRowPatFixOrd, bestRowTypFixOrd);
                if (bestErrFixOrd < bestErr) {
                    bestErr = bestErrFixOrd;
                    bestFreeSet.assign(freeSets[i].begin(), freeSets[i].end());
                    bestBoundSet.assign(boundSets[i].begin(), boundSets[i].end());
                    bestFreeSetIds[tabId] = i;
                    bestRowPat.assign(bestRowPatFixOrd.begin(), bestRowPatFixOrd.end());
                    bestRowTyp.assign(bestRowTypFixOrd.begin(), bestRowTypFixOrd.end());
                }
            }
#ifdef DEBUG
            std::cout << "best free set " << Vect2String(bestFreeSet) << std::endl;
            std::cout << "best bound set " << Vect2String(bestBoundSet) << std::endl;
            std::cout << "best row pattern " << Vect2String(bestRowPat) << std::endl;
            std::cout << "best row type " << Vect2String(bestRowTyp) << std::endl;
            std::cout << "best error " << bestErr << std::endl;
#endif
            // join approximate truth tables
            JoinAppTruTables(bestRowPat, bestRowTyp, tabId);

            // update free/bound set table
            GenFreeBoundTable(bestRowPat, bestRowTyp, tabId);
        }
    }

    void IncApprDec(int startBit, int nTry, int nComb, bool useIlp) {
        // enumerate free sets and bound sets
        std::vector < std::vector <int> > freeSets;
        std::vector < std::vector <int> > boundSets;
        if (nIn <= 10)
            EnumVarOrd(freeSets, boundSets);
        else
            DrawVarOrd(nComb, freeSets, boundSets);
        if (nComb > 0 && nComb < static_cast<int>(freeSets.size())) {
            freeSets.resize(nComb);
            boundSets.resize(nComb);
        }
#ifdef DEBUG
        assert(startBit >= 0 && startBit <= nOut - 1);
        assert(nBound + nFree == nIn);
        assert(nOut == static_cast <int> (truthTables.size()));
        assert(freeSets.size() < (1 << 20));
        assert(freeSets.size() == boundSets.size());
        std::cout << "startBit = " << startBit << std::endl;
        int accBit = (1 << nIn) * nOut;
        int appBit = ((1 << nBound) + (1 << (nFree + 1))) * (startBit + 1) + (1 << nIn) * (nOut - startBit - 1);
        std::cout << "bits of accurate lut = " << accBit << std::endl;
        std::cout << "bits of approximate lut = " << appBit << std::endl;
        std::cout << "bit saving = " << 1 - appBit / static_cast <double> (accBit) << std::endl;
        std::cout << "size of partitions = " << freeSets.size() << std::endl;
#endif
        // decompose each truth table
        appTable.assign(table.begin(), table.end());
        bestFreeSets.resize(nOut);
        bestBoundSets.resize(nOut);
        bestFreeSetIds.resize(nOut);
        freeTable.resize(1 << (nFree + 1), 0);
        boundTable.resize(1 << nBound, 0);
        for (int trial = 0; trial < nTry; ++trial) {
            for (int tabId = startBit; tabId >= 0; --tabId) {
#ifdef DEBUG
                std::cout << "**************** approximating " << tabId << "-th table ****************" << std::endl;
#endif
                double bestErr = DBL_MAX;
                std::vector <bool> bestRowPat;
                std::vector <int> bestRowTyp;
                std::vector <int> & bestFreeSet = bestFreeSets[tabId];
                std::vector <int> & bestBoundSet = bestBoundSets[tabId];
                // for each free/bound set pair, find the best decomposition
                boost::timer::progress_display procBar(freeSets.size());
                for (int i = 0; i < static_cast <int> (freeSets.size()); ++i) {
                    std::vector <bool> bestRowPatFixOrd(1 << nBound);
                    std::vector <int> bestRowTypFixOrd((1 << nFree));
                    double bestErrFixOrd = useIlp
                        ? SolveIncOpt(tabId, freeSets[i], boundSets[i], bestRowPatFixOrd, bestRowTypFixOrd)
                        : SolveIncOptADMM(tabId, freeSets[i], boundSets[i], bestRowPatFixOrd, bestRowTypFixOrd);
                    if (bestErrFixOrd < bestErr) {
                        bestErr = bestErrFixOrd;
                        bestFreeSet.assign(freeSets[i].begin(), freeSets[i].end());
                        bestBoundSet.assign(boundSets[i].begin(), boundSets[i].end());
                        bestFreeSetIds[tabId] = i;
                        bestRowPat.assign(bestRowPatFixOrd.begin(), bestRowPatFixOrd.end());
                        bestRowTyp.assign(bestRowTypFixOrd.begin(), bestRowTypFixOrd.end());
                    }
                    ++procBar;
                }
#ifdef DEBUG
                std::cout << "best free set " << Vect2String(bestFreeSet) << std::endl;
                std::cout << "best bound set " << Vect2String(bestBoundSet) << std::endl;
                std::cout << "best row pattern " << Vect2String(bestRowPat) << std::endl;
                std::cout << "best row type " << Vect2String(bestRowTyp) << std::endl;
                std::cout << "best error " << bestErr << std::endl;
#endif
                // join approximate truth tables
                JoinAppTruTables(bestRowPat, bestRowTyp, tabId);
                // update free/bound set table
                GenFreeBoundTable(bestRowPat, bestRowTyp, tabId);
            }
            std::cout << "trial = " << trial << std::endl;
            Eval();
        }
    }

    // void NonDisjIncApprDec(int startBit) {
    //     // enumerate free sets and bound sets
    //     std::vector < std::vector <int> > freeSets;
    //     std::vector < std::vector <int> > boundSets;
    //     NonDisjEnumVarOrd(freeSets, boundSets);
    //     assert(startBit >= 0 && startBit <= nOut - 1);
    //     assert(nBound + nFree - 1 == nIn);
    //     assert(nOut == static_cast <int> (truthTables.size()));
    //     assert(freeSets.size() < (1 << 20));
    //     assert(freeSets.size() == boundSets.size());
    //     std::cout << "startBit = " << startBit << std::endl;
    //     int accBit = (1 << nIn) * nOut;
    //     int appBit = ((1 << nBound) + (1 << (nFree + 1))) * (startBit + 1) + (1 << nIn) * (nOut - startBit - 1);
    //     std::cout << "bits of accurate lut = " << accBit << std::endl;
    //     std::cout << "bits of approximate lut = " << appBit << std::endl;
    //     std::cout << "bit saving = " << 1 - appBit / static_cast <double> (accBit) << std::endl;
    //     std::cout << "size of partitions = " << freeSets.size() << std::endl;
    //     // decompose each truth table
    //     appTable.assign(table.begin(), table.end());
    //     bestFreeSets.resize(startBit + 1);
    //     bestBoundSets.resize(startBit + 1);
    //     for (int tabId = startBit; tabId >= 0; --tabId) {
    //         std::cout << "**************** approximating " << tabId << "-th table ****************" << std::endl;
    //         double bestErr = DBL_MAX;
    //         std::vector <bool> bestRowPat;
    //         std::vector <int> bestRowTyp;
    //         std::vector <int> & bestFreeSet = bestFreeSets[tabId];
    //         std::vector <int> & bestBoundSet = bestBoundSets[tabId];
    //         // for each free/bound set pair, find the best decomposition
    //         boost::timer::progress_display procBar(freeSets.size());
    //         for (int i = 0; i < static_cast <int> (freeSets.size()); ++i) {
    //             std::vector <bool> bestRowPatFixOrd(1 << nBound);
    //             std::vector <int> bestRowTypFixOrd((1 << nFree));
    //             double bestErrFixOrd = NonDisjSolveIncOpt(tabId, freeSets[i], boundSets[i], bestRowPatFixOrd, bestRowTypFixOrd);
    //             if (bestErrFixOrd < bestErr) {
    //                 bestErr = bestErrFixOrd;
    //                 bestFreeSet.assign(freeSets[i].begin(), freeSets[i].end());
    //                 bestBoundSet.assign(boundSets[i].begin(), boundSets[i].end());
    //                 bestFreeSetIds[tabId] = i;
    //                 bestRowPat.assign(bestRowPatFixOrd.begin(), bestRowPatFixOrd.end());
    //                 bestRowTyp.assign(bestRowTypFixOrd.begin(), bestRowTypFixOrd.end());
    //             }
    //             ++procBar;
    //         }
    //         JoinAppTruTables(bestRowPat, bestRowTyp, tabId);
    //         std::cout << "best free set " << Vect2String(bestFreeSet) << std::endl;
    //         std::cout << "best bound set " << Vect2String(bestBoundSet) << std::endl;
    //         std::cout << "best row pattern " << Vect2String(bestRowPat) << std::endl;
    //         std::cout << "best row type " << Vect2String(bestRowTyp) << std::endl;
    //         std::cout << "best error " << bestErr << std::endl;
    //     }
    // }

    void FineGrainedApprDec() {
        assert(nBound + nFree == nIn);
        assert(nBound <= 10 && nFree <= 10);

        int accBit = (1 << nIn) * nOut;
        int startBit = nOut - 1;
        int appBit = ((1 << nBound) + (1 << (nFree + 1))) * (startBit + 1) + (1 << nIn) * (nOut - startBit - 1);
        std::cout << "bits of accurate lut = " << accBit << std::endl;
        std::cout << "bits of approximate lut = " << appBit << std::endl;
        std::cout << "bit saving = " << 1 - appBit / static_cast <double> (accBit) << std::endl;

        // enumerate free sets and bound sets
        std::vector < std::vector <int> > freeSets;
        std::vector < std::vector <int> > boundSets;
        EnumVarOrd(freeSets, boundSets);

        // decompose value table
        assert(nOut == static_cast <int> (truthTables.size()));
        bestFreeSets.resize(nOut);
        bestBoundSets.resize(nOut);
        appTable.resize(table.size(), 0);
        freeTable.resize(1 << (nFree + 1), 0);
        boundTable.resize(1 << nBound, 0);

        // for each free/bound set pair, find the best decomposition
        double bestErr = DBL_MAX;
        std::vector <bool> bestRowPat;
        std::vector <int> bestRowTyp;
        assert(freeSets.size() < static_cast <int> (INT_MAX));
        boost::timer::progress_display procBar(freeSets.size());
        time_t st = time(NULL);
        for (int i = 0; i < static_cast <int> (freeSets.size()); ++i) {
            std::vector <bool> bestRowPatFixOrd(nOut << nBound);
            std::vector <int> bestRowTypFixOrd((nOut << nFree));
            double bestErrFixOrd = SolveFineGrained(freeSets[i], boundSets[i], bestRowPatFixOrd, bestRowTypFixOrd);
            if (bestErrFixOrd < bestErr) {
                bestErr = bestErrFixOrd;
                bestFreeSets[0].assign(freeSets[i].begin(), freeSets[i].end());
                bestBoundSets[0].assign(boundSets[i].begin(), boundSets[i].end());
                bestRowPat.assign(bestRowPatFixOrd.begin(), bestRowPatFixOrd.end());
                bestRowTyp.assign(bestRowTypFixOrd.begin(), bestRowTypFixOrd.end());
            }
            ++procBar;
            std::cout << "current best error = " << bestErr << std::endl;
            time_t ed = time(NULL);
            std::cout << "runtime = " << ed - st << "s" << std::endl;
            break;
        }

        std::cout << "best free set " << Vect2String(bestFreeSets[0], 0, ",") << std::endl;
        std::cout << "best bound set " << Vect2String(bestBoundSets[0], 0, ",") << std::endl;
        std::cout << "best row pattern " << Vect2String(bestRowPat, 0, ",") << std::endl;
        std::cout << "best row type " << Vect2String(bestRowTyp, 0, ",") << std::endl;
        std::cout << "best error " << bestErr << std::endl;

        // patch free/bound set
        for (int i = 1; i < nOut; ++i) {
            bestFreeSets[i].assign(bestFreeSets[0].begin(), bestFreeSets[0].end());
            bestBoundSets[i].assign(bestBoundSets[0].begin(), bestBoundSets[0].end());
        }

        // // generate approximate table
        // GenAppTableAllBitsFine(bestRowPat, bestRowTyp);

        // // generate free/bound table
        // GenFreeBoundTableAllBitsFine(bestRowPat, bestRowTyp);
    }

    void GetBoundFromFree(std::vector < std::vector <int> > & freeSets, std::vector < std::vector <int> > & boundSets) {
        boundSets.resize(freeSets.size());
        for (int k = 0; k < static_cast <int> (freeSets.size()); ++k) {
            auto & freeSet = freeSets[k];
            auto & boundSet = boundSets[k];
            assert(static_cast <int> (freeSet.size()) == nFree);
            int i = 0, j = 0;
            while (i < nIn) {
                if (j < nFree) {
                    if (i != freeSet[j])
                        boundSet.emplace_back(i);
                    else
                        ++j;
                }
                else
                    boundSet.emplace_back(i);
                ++i;
            }
        }
    }

    void GetBoundFromFree(const std::vector <int> & freeSet, std::vector <int> & boundSet, int overlId) {
        assert(static_cast <int> (freeSet.size()) == nFree);
        assert(overlId < nFree);
        boundSet.clear();
        boundSet.reserve(nBound);
        int i = 0, j = 0;
        while (i < nIn) {
            if (i == freeSet[overlId]) {
                boundSet.emplace_back(i);
            }
            if (j < nFree) {
                if (i != freeSet[j])
                    boundSet.emplace_back(i);
                else
                    ++j;
            }
            else
                boundSet.emplace_back(i);
            ++i;
        }
    }

    void FineGrainedApprDecPro(int startBit) {
#ifdef DEBUG
        assert(nBound + nFree == nIn);
        assert(nBound <= 10 && nFree <= 10);
        assert(nOut == static_cast <int> (truthTables.size()));
        int accBit = (1 << nIn) * nOut;
        int appBit = ((1 << nBound) + (1 << (nFree + 1))) * (startBit + 1) + (1 << nIn) * (nOut - startBit - 1);
        std::cout << "bits of accurate lut = " << accBit << std::endl;
        std::cout << "bits of approximate lut = " << appBit << std::endl;
        std::cout << "bit saving = " << 1 - appBit / static_cast <double> (accBit) << std::endl;
#endif
        // enumerate free sets and bound sets
        std::vector < std::vector <int> > freeSets = {
            {1, 2, 3}, {1, 2, 4}, {0, 1, 4}, {0, 3, 4}, {1, 2, 4}, {1, 4, 5}, {2, 3, 5}, {0, 1, 2}
        };
        // std::vector < std::vector <int> > freeSets = {
        //     {0, 1, 5},  {0, 1, 5},  {0, 3, 4},  {0, 4, 5},  {1, 2, 3},  {0, 4, 5},  {0, 2, 5}
        // };
        // std::vector < std::vector <int> > freeSets = {
        //     {1, 4, 5}, {1, 2, 4}, {2, 3, 4}, {0, 3, 5}, {2, 3, 5}, {0, 2, 5}
        // };
        // std::vector < std::vector <int> > freeSets = {
        //     {1, 2, 4}, {0, 1, 2}, {0, 2, 5}, {0, 2, 3}, {0, 1, 4}
        // };
        // std::vector < std::vector <int> > freeSets = {
        //     {0, 3, 4}, {0, 2, 3}, {1, 4, 5}, {1, 4, 5}
        // };
        // std::vector < std::vector <int> > freeSets = {
        //     {0, 2, 3}, {0, 1, 5}, {0, 1, 3}
        // };
        // std::vector < std::vector <int> > freeSets = {
        //     {1, 2, 5}, {0, 1, 4}
        // };
        // std::vector < std::vector <int> > freeSets = {
        //     {1, 2, 5}
        // };
        std::vector < std::vector <int> > boundSets;
        GetBoundFromFree(freeSets, boundSets);

        for (int k = 0; k < static_cast <int> (freeSets.size()); ++k)
            std::cout << Vect2String(freeSets[k], 0, ",") << "," << Vect2String(boundSets[k], 0, ",") << std::endl;

        // find the best decomposition
        std::vector <bool> bestRowPat((startBit + 1) << nBound);
        std::vector <int> bestRowTyp((startBit + 1) << nFree);
        SolveFineGrainedPro(startBit, freeSets, boundSets, bestRowPat, bestRowTyp);
    }

    void Eval() {
#ifdef DEBUG
        assert(table.size() == appTable.size());
#endif
        double med = 0.0;
        for (int i = 0; i < static_cast <int> (table.size()); ++i) {
#ifdef DEBUG
            assert(table[i] <= (static_cast <uint64_t> (1) << 62));
            assert(appTable[i] <= (static_cast <uint64_t> (1) << 62));
#endif
            med += abs(static_cast <int64_t> (table[i]) - static_cast <int64_t> (appTable[i])) * prob[i];
        }
        // med /= maxOutVal;
        std::cout << "med = " << med << std::endl;
        std::cout << "nmed = " << med / static_cast <double> (maxOutVal) * 100 << "%" << std::endl;

        double mse = 0.0;
        for (int i = 0; i < static_cast <int> (table.size()); ++i) {
#ifdef DEBUG
            assert(table[i] <= (static_cast <uint64_t> (1) << 62));
            assert(appTable[i] <= (static_cast <uint64_t> (1) << 62));
#endif
            int64_t diff = static_cast <int64_t> (table[i]) - static_cast <int64_t> (appTable[i]);
            // std::cout << i << "," << table[i] << "," << appTable[i] << std::endl;
            mse +=  prob[i] * diff * diff;
        }
        const double epsilon = 1e-30;
        double psnr = 10 * log10(pow(static_cast <double> (maxOutVal), 2) / (mse + epsilon));
        std::cout << "mse = " << mse << std::endl;
        std::cout << "psnr = " << psnr << std::endl;
    }

    void Output(std::string output, std::string fileName) {
        // prepare output path
        if (output[output.length() - 1] != '/')
            output += "/";
        std::filesystem::path outputPath(output);
        if (!std::filesystem::exists(outputPath))
            std::filesystem::create_directories(outputPath);
        assert(std::filesystem::exists(outputPath));

        // dump
        DumpHardwPar(output + fileName);

        // output approximate table
        // std::string outFile = output + fileName;
        // FILE * file = fopen(outFile.c_str(), "w");
        // assert(file != nullptr);
        // for (int i = 0; i < nAddr; ++i)
        //     fprintf(file, "%d %lu\n", i, appTable[i]);
        // fclose(file);
    }

    void DumpHardwPar(std::string outPrefix) const {
        FILE * fp = nullptr;
        fp = fopen((outPrefix + "_accTab.txt").c_str(), "w");
        assert(fp != nullptr);
        for (int i = 0; i < static_cast <int> (table.size()); ++i)
            fprintf(fp, "%d %ld\n", i, table[i]);
        fclose(fp);

        fp = fopen((outPrefix + "_appTab.txt").c_str(), "w");
        assert(fp != nullptr);
        for (int i = 0; i < static_cast <int> (appTable.size()); ++i)
            fprintf(fp, "%d %ld\n", i, appTable[i]);
        fclose(fp);

        assert(nIn <= 20); // nIn * log2(nIn) <= 127
        assert(nOut == static_cast <int> (bestFreeSets.size()));
        assert(nOut == static_cast <int> (bestBoundSets.size()));
        int nSel = ceil(log2(nIn));
        fp = fopen((outPrefix + "_appSel.txt").c_str(), "w");
        typedef boost::multiprecision::uint128_t bigInt;
        for (int i = 0; i < nOut; ++i) {
            if (bestFreeSets[i].size() == 0)
                break;
            bigInt selSig = 0;
            for (int j = nFree - 1; j >= 0; --j) {
                selSig <<= nSel;
                selSig |= bestFreeSets[i][j];
            }
            for (int j = nBound - 1; j >= 0; --j) {
                selSig <<= nSel;
                selSig |= bestBoundSets[i][j];
            }
            uint64_t selSigH = static_cast <uint64_t> (selSig >> 64);
            uint64_t selSigL = static_cast <uint64_t> (selSig & UINT64_MAX);
            if (selSigH)
                fprintf(fp, "%lu%lu\n", selSigH, selSigL);
            else
                fprintf(fp, "%lu\n", selSigL);
        }
        fclose(fp);

        fp = fopen((outPrefix + "_freeTab.txt").c_str(), "w");
        assert(fp != nullptr);
        for (int i = 0; i < static_cast <int> (freeTable.size()); ++i)
            fprintf(fp, "%d %ld\n", i, freeTable[i]);
        fclose(fp);

        fp = fopen((outPrefix + "_boundTab.txt").c_str(), "w");
        assert(fp != nullptr);
        for (int i = 0; i < static_cast <int> (boundTable.size()); ++i)
            fprintf(fp, "%d %ld\n", i, boundTable[i]);
        fclose(fp);
    }

    void EvalRound()
    {
        assert(nOut <= 20);
        int accBit = (1 << nIn) * nOut;
        for (int nRound = 1; nRound < nOut; ++nRound) {
            int appBit = (1 << nIn) * (nOut - nRound);
            std::cout << "round " << nRound << " bits" << std::endl;
            std::cout << "bit saving = " << 1 - appBit / static_cast <double> (accBit) << std::endl;
            double med = 0;
            std::vector <uint64_t> newVal(table.size());
            for (int i = 0; i < static_cast <int> (table.size()); ++i) {
                newVal[i] = table[i] >> nRound;
                int carry = GetBit(table[i], nRound - 1);
                newVal[i] += carry;
                newVal[i] <<= nRound;
                if (static_cast <int> (newVal[i]) > maxOutVal)
                    newVal[i] = maxOutVal;
                med += abs(static_cast <int64_t> (table[i]) - static_cast <int64_t> (newVal[i])) * prob[i];
            }
            // std::cout << Vect2String(table, 0, ",") << std::endl;
            // std::cout << Vect2String(newVal, 0, ",") << std::endl;
            std::cout << "med = " << med << std::endl;
            std::cout << "nmed = " << med / static_cast <double> (maxOutVal) * 100 << "%" << std::endl;
            double mse = 0.0;
            for (int i = 0; i < static_cast <int> (table.size()); ++i) {
                int64_t diff = static_cast <int64_t> (table[i]) - static_cast <int64_t> (newVal[i]);
                mse +=  prob[i] * diff * diff;
            }
            const double epsilon = 1e-30;
            double psnr = 10 * log10(pow(static_cast <double> (maxOutVal), 2) / (mse + epsilon));
            std::cout << "mse = " << mse << std::endl;
            std::cout << "psnr = " << psnr << std::endl;
        }
    }

};


#endif
