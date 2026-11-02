#ifndef BIT_OPER_HPP
#define BIT_OPER_HPP


#include "headers.h"


const uint64_t TABLE_LEN = (1 << 22);


class Bit_Count_t
{
private:
    int table[TABLE_LEN];

    Bit_Count_t & operator == (const Bit_Count_t &);
    Bit_Count_t(Bit_Count_t &);
    Bit_Count_t(Bit_Count_t &&);

public:
    explicit Bit_Count_t (void) {
        std::cout << "constructor" << std::endl;
        for (uint64_t i = 0; i < TABLE_LEN; ++i) {
            table[i] = 0;
            int x = i;
            while (x) {
                x = x & (x - 1);
                ++table[i];
            }
        }
    }


    ~Bit_Count_t() {
    }


    inline int GetOneNum(uint64_t n) {
        return table[n >> 44] + table[(n >> 22) & 0x3fffff] + table[n & 0x3fffff];
    }
};


template <typename T>
concept INT_NUMB = (
        std::is_same_v<T, uint8_t> ||
        std::is_same_v<T, uint16_t> ||
        std::is_same_v<T, uint32_t> ||
        std::is_same_v<T, uint64_t> ||
        std::is_same_v<T, int32_t>
);


template <INT_NUMB T>
inline void SetBit(T & x, uint32_t bit) {
    x |= (static_cast <T> (1) << static_cast <T> (bit));
}


template <INT_NUMB T>
inline void ResetBit(T & x, uint32_t bit) {
    x &= ~(static_cast <T> (1) << static_cast <T> (bit));
}


template <INT_NUMB T>
inline bool GetBit(const T x, uint32_t bit) {
    return static_cast <bool> ((static_cast <T> (x) >> static_cast <T> (bit)) & static_cast <T> (1));
}


template <INT_NUMB T>
inline int CountOneNumb(const T x) {
    static Bit_Count_t bitCnt;
    return static_cast <int> (bitCnt.GetOneNum(x));
}


template <INT_NUMB T>
inline std::string Numb2String(const T x, int msbPos) {
    assert(msbPos < 64);
    std::ostringstream oss("");
    std::vector <int> bits(64, 0);
    int count = 0;
    uint64_t numb = static_cast <uint64_t> (x);
    while (numb) {
        bits[count] = numb & static_cast <uint64_t> (1);
        numb >>= 1;
        ++count;
    }
    for (int i = msbPos; i >= 0; --i)
        oss << bits[i];
    return oss.str();
}


#endif
