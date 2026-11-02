#ifndef UTILS_HPP
#define UTILS_HPP


#include "headers.h"


inline int ExecSystComm(const std::string cmd)
{
    pid_t status;
    std::cout << "Execute system command: " << cmd << std::endl;
    status = system(cmd.c_str());
    if (status == -1)
    {
        std::cout << "System error!" << std::endl;
        return 0;
    }
    else
    {
        std::cout << "Exit status value = " << status << std::endl;
        if (WIFEXITED(status))
        {
            if (!WEXITSTATUS(status))
            {
                std::cout << "Run shell script successfully." << std::endl;
                return 1;
            }
            else
            {
                std::cout << "Run shell script fail, script exit code: " << WEXITSTATUS(status) << std::endl;
                return 0;
            }
        }
        else
        {
            std::cout << "Exit status = " << WEXITSTATUS(status) << std::endl;
            return 0;
        }
    }
}


template <typename T>
concept NUMB = (
        std::is_same_v <T, bool> ||
        std::is_same_v <T, int> ||
        std::is_same_v <T, double> ||
        std::is_same_v <T, uint64_t>
);


template <NUMB T>
inline std::string Vect2String(const std::vector <T> & vect, bool isRev = 0, std::string sep = ",") {
    std::ostringstream oss("");
    oss << "[";
    if (isRev) {
        for (int i = vect.size() - 1; i >= 0; --i) {
            if (i != 0)
                oss << vect[i] << sep;
            else
                oss << vect[i];
        }
    }
    else {
        for (int i = 0; i < static_cast <int> (vect.size()); ++i) {
            if (i != static_cast <int> (vect.size()) - 1)
                oss << vect[i] << sep;
            else
                oss << vect[i];
        }
    }
    oss << "]";
    return oss.str();
}


template <NUMB T>
inline std::string VectSqErr2String(const std::vector <T> & vect0, const std::vector <T> & vect1, bool isRev, std::string sep) {
    assert(vect0.size() == vect1.size());
    std::ostringstream oss("");
    if (isRev) {
        for (int i = vect0.size() - 1; i >= 0; --i) {
            if (i != 0)
                oss << (vect0[i] - vect1[i]) *  (vect0[i] - vect1[i]) << sep;
            else
                oss << (vect0[i] - vect1[i]) *  (vect0[i] - vect1[i]);
        }
    }
    else {
        for (int i = 0; i < vect0.size(); ++i) {
            if (i != vect0.size() - 1)
                oss << (vect0[i] - vect1[i]) *  (vect0[i] - vect1[i]) << sep;
            else
                oss << (vect0[i] - vect1[i]) *  (vect0[i] - vect1[i]);
        }
    }
    return oss.str();
}


template <NUMB T>
inline std::string List2String(const std::list <T> & lis, std::string sep) {
    std::ostringstream oss("");
    for (auto elem: lis)
        oss << elem << sep;
    return oss.str();
}


#endif
