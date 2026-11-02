#include "headers.h"
#include "cmdline.hpp"
#include "utils.hpp"
#include "lut.hpp"


cmdline::parser CommPars(int argc, char * argv[])
{
    cmdline::parser option;
    option.add <std::string> ("input", '\0', "input file", false, "./benchmark/images/lena_64x64_gray.txt");
    option.add <std::string> ("output", '\0', "output path", false, "./apprBenchmark/");
    option.add <int> ("nFree", '\0', "#free variables", false, 3);
    // option.add <int> ("nBound", '\0', "#bound variables", false, 3);
    option.add <int> ("startBit", '\0', "start bit for incremental approximate decomposition", false, 3);
    option.add <int> ("nTry", '\0', "trial time", false, 5);
    option.add <int> ("nComb", '\0', "#maximum combinations of free/bound sets", false, 1000);
    option.add <int> ("admmReInitTime", '\0', "ADMM re-initialize time", false, 30);
    option.add <int> ("admmStuckTol", '\0', "ADMM stucking tolerance", false, 3);
    option.add <unsigned> ("seed", '\0', "seed", false, 0);
    option.add <double> ("MIPTime", '\0', "time limit for one Gurobi model", false, 3600.0);
    option.add <double> ("MIPGap", '\0', "relative Gurobi MIP gap", false, 0.1);
    option.add <std::string> ("mode", '\0', "mode: decompose, round", false, "decompose");
    option.add <std::string> ("solver", '\0', "solver: separate-ilp, separate-legacy, joint-ilp, joint-heuristic", false, "joint-heuristic");
    option.parse_check(argc, argv);
    return option;
}


int main(int argc, char * argv[])
{
    // command line parser
    cmdline::parser option = CommPars(argc, argv);
    std::string input = option.get <std::string> ("input");
    std::string output = option.get <std::string> ("output");
    int nFree = option.get <int> ("nFree");
    // int nBound = option.get <int> ("nBound");
    int startBit = option.get <int> ("startBit");
    int nTry = option.get <int> ("nTry");
    int nComb = option.get <int> ("nComb");
    int admmReInitTime = option.get <int> ("admmReInitTime");
    int admmStuckTol = option.get <int> ("admmStuckTol");
    unsigned seed = option.get <unsigned> ("seed");
    double MIPTime = option.get <double> ("MIPTime");
    double MIPGap = option.get <double> ("MIPGap");
    std::string mode = option.get <std::string> ("mode");
    std::string solver = option.get <std::string> ("solver");

    // approximate decomposition
    std::cout << "mode = " << mode << std::endl;
    std::cout << "solver = " << solver << std::endl;
    if (seed == 0) {
        srand(time(nullptr));
        seed = static_cast <unsigned> (rand());
    }
    std::cout << "seed = " << seed << std::endl;
    std::filesystem::path inputPath(input);
    assert(exists(inputPath));
    LUT_t lut(input, nFree, admmReInitTime, admmStuckTol, seed, MIPTime, MIPGap);
    clock_t st = clock();
    if (mode == "decompose") {
        lut.ApprDec(startBit, nTry, nComb, solver);
        lut.Eval();
        std::ostringstream oss;
        oss << "input_" << inputPath.stem() << "_nFree_" << nFree << "_startBit_" << startBit << "_nTry_" << nTry << "_nComb_" << nComb << "_aR_" << admmReInitTime << "_aS_" << admmStuckTol << "_seed_" << lut.GetSeed();
        // std::cout << oss.str() << std::endl;
        lut.Output(output, oss.str());
    }
    else
        lut.EvalRound();
    std::cout << "time = " << static_cast<double>(clock() - st) / CLOCKS_PER_SEC << " s" << std::endl;

    return 0;
}
