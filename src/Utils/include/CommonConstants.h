#ifndef COMMONCONSTANTS_H
#define COMMONCONSTANTS_H

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>
// #include <iostream>
#ifndef PACK_STRUCT
    #define PACK_STRUCT  __attribute__((packed))
#endif

#define CONST_C0 (299792458.) //       % speed of light in the atmosphere (m/s)
// #define ELEMENT3D(AA,BB,CC,MAT) *( reinterpret_cast<double*>((MAT).data) + ((AA)*(MAT).step[0] / sizeof(double)) + ((BB)*(MAT).step[1] / sizeof(double)) + (CC) )
// #define ELEMENT1D(__IDX__,MAT) *(reinterpret_cast<double*>((MAT).data) + (__IDX__))
// #define ELEMENT2D(AA,BB,MAT) *(reinterpret_cast<double*>((MAT).data) + ((AA)*(MAT).size[1]) + (BB) )
// #define ELEMENT1D_CPL(AA,MAT) *(reinterpret_cast< std::complex<double>* >((MAT).data) + (AA) )
// #define ELEMENT2D_CPL(AA,BB,MAT) (*( reinterpret_cast< std::complex<double>* >((MAT).data) + ((AA)*(MAT).size[1]) + (BB) ))
// #define ELEMENT2D_F(AA,BB,MAT) *(reinterpret_cast<float*>((MAT).data) + ((AA)*(MAT).size[1]) + (BB) )
// #define ELEMENT2D_U8(AA,BB,MAT) *(reinterpret_cast<uint8_t*>((MAT).data) + ((AA)*(MAT).size[1]) + (BB) )
// #define ELEMENT2D_U8C4(AA,BB,MAT) *(reinterpret_cast<uint32_t*>((MAT).data) + ((AA)*(MAT).size[1]) + (BB) )


#define ELEMENT3D(AA,BB,CC,MAT) (MAT[(AA) + (BB) * dim0 + (CC) * dim0 * dim1])
#define ELEMENT1D(__IDX__,MAT) (MAT.data()[__IDX__])
#define ELEMENT2D(AA,BB,MAT) (MAT.data()[(AA) + (BB) * (MAT).rows()])
#define ELEMENT1D_CPL(AA,MAT) (MAT.data()[AA])
#define ELEMENT2D_CPL(AA,BB,MAT) (MAT.data()[(AA) + (BB) * (MAT).rows()])
#define ELEMENT2D_F(AA,BB,MAT) (MAT.data()[(AA) + (BB) * (MAT).rows()])
#define ELEMENT2D_U8(AA,BB,MAT) (MAT.data()[(AA) + (BB) * (MAT).rows()])
#define ELEMENT2D_U8C4(AA,BB,MAT) (&(MAT.data()[(AA) * 4 + (BB) * (MAT).rows()]))

#define ELEMENT3D_FAST(AA,BB,CC,DATA,DIM0,DIM1) (DATA[(AA) + (BB) * (DIM0) + (CC) * (DIM0) * (DIM1)])
#define ELEMENT1D_FAST(__IDX__,DATA) (DATA[__IDX__])
#define ELEMENT2D_FAST(AA,BB,DATA,ROWS) (DATA[(AA) + (BB) * (ROWS)])
#define ELEMENT1D_CPL_FAST(AA,DATA) (DATA[AA])
#define ELEMENT2D_CPL_FAST(AA,BB,DATA,ROWS) (DATA[(AA) + (BB) * (ROWS)])
#define ELEMENT2D_F_FAST(AA,BB,DATA,ROWS) (DATA[(AA) + (BB) * (ROWS)])
#define ELEMENT2D_U8_FAST(AA,BB,DATA,ROWS) (DATA[(AA) + (BB) * (ROWS)])
#define ELEMENT2D_U8C4_FAST(AA,BB,DATA,ROWS) (&(DATA[(AA) * 4 + (BB) * (ROWS)]))



#define CONST_ACCEPTABLE_MAX_DIFF_ERROR 1.e-10
#define CONST_ACCEPTABLE_SUM_DIFF_ERROR 1.e-5


template<typename T>
std::vector<T>& operator<<(std::vector<T>& vec, const T& value) {
    vec.push_back(value);
    return vec;
}

template<typename T>
bool removeOne(std::vector<T>& _vec, const T& _obj)
{
    const auto _iter = std::find(_vec.begin(), _vec.end(), _obj);
    if(_vec.end() != _iter)
    {
        _vec.erase(_iter);
        return true;
    }
    return false;
}
template<typename T>
bool contains(const std::vector<T>& _vec, const T& _obj)
{
    return (std::find(_vec.begin(), _vec.end(), _obj) != _vec.end());
}

template<typename T1, typename T2>
std::vector<T1> map_keys(const std::unordered_map<T1, T2>& _hash)
{
    std::vector<T1> _list;
    for(const auto& _iter: _hash)
    {
        _list.push_back(_iter.first);
    }
    return _list;
}

template<typename T1, typename T2>
std::vector<T2> map_values(const std::unordered_map<T1, T2>& _hash)
{
    std::vector<T1> _list;
    for(const auto& _iter: _hash)
    {
        _list.push_back(_iter.second);
    }
    return _list;
}

/**
 *
 * converts input to string
 */
#define STRING(s) #s

static void assert_condition(bool condition, const char* condition_str, const std::string& message)
{
    if(!condition) {
        // std::cerr << "Assertion failed: " << condition_str << " | Message: " << message << std::endl;
        throw std::runtime_error("Assertion failed: " + std::string(condition_str) + " | " + message);
    }
}
#define assertm(COND, MSG) assert_condition((COND), #COND, MSG)

class CommonConstants
{
public:
    CommonConstants() {}
    virtual ~CommonConstants() {}

private:

public:

    
};




#endif // COMMONCONSTANTS_H
