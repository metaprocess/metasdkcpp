#ifndef __DEFINITIONS_H__
#define __DEFINITIONS_H__
#include <math.h>
#include <iomanip>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <mutex>
// #include "structures.h"
#include <stdint.h>
#include <memory>
#include <sstream>
#include "CommonConstants.h"



struct uint24_t {
    uint32_t v:24;
    uint24_t(){};
    uint24_t(const uint32_t& _v): v(_v){}
    bool operator==(const uint24_t& inp) const { return v == inp.v; }
    bool operator>(const uint24_t& inp) const { return v > inp.v; }
    bool operator<(const uint24_t& inp) const { return v < inp.v; }
    bool operator>=(const uint24_t& inp) const { return v >= inp.v; }
    bool operator<=(const uint24_t& inp) const { return v <= inp.v; }
    // Prefix increment operator
    uint24_t& operator++() {
        if (v < 0xFFFFFF) { // Ensure the value does not exceed 24 bits
            ++v;
        }
        return *this;
    }
    // Postfix increment operator
    uint24_t operator++(int) {
        uint24_t temp = *this;
        ++(*this); // Use the prefix operator to increment the value
        return temp; // Return the value before increment
    }
}PACK_STRUCT;

static std::vector<std::string> split_str(const std::string& _str, const char &_sep)
{
    std::vector<std::string> _list_str;
    std::stringstream ss(_str);

    // Use while loop to check the getline() function condition.
    std::string str;
    while(std::getline(ss, str, _sep))
    {
        // `str` is used to store the token string while 'separator(_sep)' is used as the delimiter.
        _list_str.push_back(str);
    }

    return _list_str;
}


static const int g_days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
static const int j_days_in_month[12] = {31, 31, 31, 31, 31, 31, 30, 30, 30, 30, 30, 29};
static const char* j_month_names[12] = {"Farvardin", "Ordibehesht", "Khordad", "Tir", "Mordad", "Shahrivar", "Mehr", "Aban", "Azar", "Dey", "Bahman", "Esfand"};

static const struct {
    const double SemimajorAxis{6.378137e+6};
    const double SemiminorAxis{6.35675231424518e+6};
    const double SquaredEccentricity{((SemimajorAxis*SemimajorAxis) - (SemiminorAxis*SemiminorAxis)) / (SemimajorAxis*SemimajorAxis)};
    const double Flattening{(SemimajorAxis - SemiminorAxis) / SemimajorAxis};
} wgs84Meter;


//#define CONST_SERIAL_DATA_LEN 7

#ifndef nullptr
#define nullptr NULL
#endif

static const char* global_timestamp = __TIMESTAMP__;

#define DEGREE_TO_RADIAN_UNIT  (M_PI/180.0)
#define RADIAN_TO_DEGREE_UNIT  (180.0/M_PI)


struct systemInfo
{
    systemInfo()
    {
        bzero(this, sizeof(systemInfo));
    }
    // in Mega Byte
    int TotalMemory_MB;
    int AvailableMemory_MB;
    int ProcId;
    double physMemUsed_MB;
    double totalSwap;
    double freeSwap;
    double usedSwap;
    //Virtual Memory currently used by current process in Mb
    double vmProc;
    //Physical Memory currently used by current process in Mb
    double PhysMem;
    int Numthread;
    double SysCpuPercent;
    double ProcCpuPercent;
    double coreTemp;
    int64_t uptime_seconds;
    uint64_t disk_space_free;
    uint64_t disk_space_total;
    int64_t throughput_link;
};



class NullBuffer : public std::streambuf
{
public:
  int overflow(int c) { return c; }
};

#if 1 == DEV_DEBUG
    #define COND_CERR std::cerr
    #define DEBUG_FUNCTION_LINE std::cerr << Utils::getTickCountMs() << "::"<< __FILE__ << "::" << __FUNCTION__ << "::" << __LINE__ << "\n";
#else
    #define COND_CERR null_stream
    #define DEBUG_FUNCTION_LINE
#endif

#define DEBUG_FUNCTION_LINE_PERMANENT std::cerr << Utils::getTickCountMs() << "::"<< __FILE__ << "::" << __FUNCTION__ << "::" << __LINE__ << "\n";



class DataContent
{
private:
    size_t content_size{};
    uint8_t* content_data{};

public:
    explicit DataContent(const char* name = nullptr);

    DataContent(const DataContent& dc);
    DataContent(DataContent&& _other) noexcept;

    ~DataContent();
    inline size_t size() const
    {
        return content_size;
    }

    void reserve(const size_t& res_size);

    void set_data(const void* input_data, const size_t& input_size);

    void insert_to_beginning(const void* input_data, const size_t& input_size);

    void add_to_end(const void* input_data, const size_t& input_size);

    void clear();

    inline const uint8_t* begin() const
    {
        return content_data;
    }

    inline const uint8_t* end() const
    {
        return (content_data + content_size);
    }

    inline uint8_t* data() const
    {
        return content_data;
    }
    inline const uint8_t* data_const() const
    {
        return content_data;
    }

    inline uint8_t& operator[](const size_t& index)
    {
        if (index >= content_size)
        {
            throw std::out_of_range("Index out of range");
        }
        return content_data[index];
    }

    inline const uint8_t& operator[](const size_t& index) const
    {
        if (index >= content_size)
        {
            throw std::out_of_range("Index out of range");
        }
        return content_data[index];
    }
};

enum Level
{
    Level_None = 0,
    Level_All = 1,
    Level_Info = (1 << 1),
    Level_Debug = (1 << 2),
    Level_Warning = (1 << 3),
    Level_Error = (1 << 4),
    Level_Critical = (1 << 5)
};

struct ComplexType
{
    ComplexType()
    {
        bzero(this, sizeof(ComplexType));
    }
    ComplexType(const double& re, const double& im)
    {
        real = re;
        imag = im;
    }
    double real;
    double imag;
    //	void operator*=(const complexType& cpl)
    //	{

    //	}
    inline void operator-=(const ComplexType& cpl)
    {
        real -= cpl.real;
        imag -= cpl.imag;
    }

    inline void operator=(const double& val)
    {
        real = imag = val;
    }
    inline void operator+=(const ComplexType& cpl)
    {
        real += cpl.real;
        imag += cpl.imag;
    }

    inline void operator*=(const ComplexType& cpl)
    {
        double re = (real * cpl.real)-(imag*cpl.imag);
        double im = real*cpl.imag + cpl.real*imag;
        real = re;
        imag = im;
    }

    inline void operator*=(const double& val)
    {
        real *= val;
        imag *= val;
    }

    inline void operator/=(const double& divi)
    {
        real /= divi;
        imag /= divi;
    }

    inline ComplexType operator*(ComplexType cpl1/*, const complexType& cpl2*/)
    {
        ComplexType ct(this->real, this->imag);
        //    ct.real = (cpl1.real * cpl2.real)-(cpl1.imag*cpl2.imag);
        //    ct.imag = cpl1.real*cpl2.imag + cpl2.real*cpl1.imag;
        ct *= cpl1;
        return ct;
    }

    inline ComplexType operator*(const double& cpl1/*, const complexType& cpl2*/)
    {
        ComplexType ct(this->real, this->imag);
        //    ct.real = (cpl1.real * cpl2.real)-(cpl1.imag*cpl2.imag);
        //    ct.imag = cpl1.real*cpl2.imag + cpl2.real*cpl1.imag;
        ct *= cpl1;
        return ct;
    }

    inline void operator=(const ComplexType& cpl)
    {
        memcpy(this, &cpl, sizeof(ComplexType));
    }

    inline double abs()
    {
        //		return JavaMath::sqrt( (real * real) + (imag * imag) );
        return sqrt( (real * real) + (imag * imag) );
    }
};

struct PointDouble
{
    double x;
    double y;

    // PointDouble()
    // {
    //     x = 0;
    //     y = 0;
    // }

    // PointDouble(const PointDouble& _other)
    // {
    //     x = _other.x;
    //     y = _other.y;
    // }

    // PointDouble(const double& _x, const double& _y)
    // {
    //     x = _x;
    //     y = _y;
    // }
    // inline PointDouble& operator=(const PointDouble& _rhs)
    // {
    //     x = _rhs.x;
    //     y = _rhs.y;
    //     return *this;
    // }
    inline PointDouble& operator*=(double rhs)
    {
        x *= rhs;
        y *= rhs;
        return *this;
    }

    inline PointDouble operator*(double rhs)
    {
        PointDouble _d;
        _d.x = x * rhs;
        _d.y = y * rhs;
        return _d;
    }

    inline bool operator==(const double& d) const
    {
        return (0. == std::abs(d - x) && 0. == std::abs(d - y));
    }

    inline bool operator==(const PointDouble& _p) const
    {
        return (0. == std::abs(_p.x - x) && 0. == std::abs(_p.y - y));
    }

    inline bool operator!=(const PointDouble& _p) const
    {
        return !(_p == *this);
    }

    inline bool operator!=(const double& d) const
    {
        return (0. != std::abs(d - x) || 0. != std::abs(d - y));
    }

    std::string toString() const
    {
        std::ostringstream oss;
        oss << "x:" << std::setprecision(10) << x << ",y:" << y;
        return oss.str();
    }

    static PointDouble from_xy(const double& _x, const double& _y)
    {
        PointDouble _pd;
        _pd.x = _x;
        _pd.y = _y;
        return _pd;
    }

}PACK_STRUCT;


struct PointPolar3D
{
    PointPolar3D()
    {
        bzero(this, sizeof(PointPolar3D));
    }
    PointPolar3D(const double& r, const double& t, const double& z)
    {
        this->r = r;
        this->t = t;
        this->z = z;
    }
    inline PointPolar3D& operator=(const PointPolar3D& point)
    {
        memcpy(this, &point, sizeof(PointPolar3D));
        return *this;
    }
    inline bool operator==(const double& d) const
    {
        return (d == r && d == t && d == z);
    }
    inline bool operator!=(const double& d) const
    {
        return (d != r || d != t || d != z);
    }
    inline bool operator!=(const PointPolar3D& d)
    {
        return (d.r != r || d.t != t || d.z != z);
    }

    bool isValid()
    {
        bool is_valid{ true };

        is_valid &= std::isnormal(r);
        is_valid &= (t >= -180 && t <= 360);
        is_valid &= std::isfinite(z);

        return is_valid;
    }

    double r;
    double t;
    double z;
}PACK_STRUCT;

struct PointDouble3DLight;

struct PointDouble3D
{
    PointDouble3D()
    {
        bzero(this, sizeof(PointDouble3D));
    }
    PointDouble3D(const double& _x, const double& _y, const double& _z)
    {
        set(_x, _y, _z);
    }
    PointDouble3D(const PointDouble3DLight& _other)
    {
        *this = _other;
    }
    ~PointDouble3D() {}
    inline PointDouble3D& operator=(const PointDouble3D& point)
    {
        memcpy(this, &point, sizeof(PointDouble3D));
        return *this;
    }
    void set(const double& _x, const double& _y, const double& _z)
    {
        this->x = _x;
        this->y = _y;
        this->z = _z;
    }
    PointDouble3D& operator=(const PointDouble3DLight& point);
    inline bool operator==(const double& d) const
    {
        return (0. == std::abs(d - x) && 0. == std::abs(d - y) && 0. == std::abs(d - z));
    }
    inline bool operator!=(const double& d) const
    {
        return (0. != std::abs(d - x) || 0. != std::abs(d - y) || 0. != std::abs(d - z));
    }
    inline bool operator!=(const PointDouble3D& d)
    {
        return (0. != std::abs(d.x - x) || 0. != std::abs(d.y - y) || 0. != std::abs(d.z - z));
    }
    friend std::ostream &operator<<( std::ostream &output, const PointDouble3D &_p )
    {
      output << std::setprecision(6) << std::fixed << "X: " << _p.x << ", Y: " << _p.y << ", Z: " << _p.z;
      return output;
    }


    inline double abs3()
    {
        return sqrt(x*x+y*y+z*z);
    }

    inline void normalize()
    {
        double absVal = abs3();
        x /= absVal;
        y /= absVal;
        z /= absVal;
    }

    const PointDouble to_point_double() const
    {
        return PointDouble::from_xy(x, y);
    }

    double y;
    double x;
    double z;
}PACK_STRUCT;

bool operator==(const double& v, const PointDouble3D& pp);

bool operator!=(const double& v, const PointDouble3D& pp);

struct PointDouble3DLight
{
    PointDouble3DLight()
    {
        bzero(this, sizeof(PointDouble3DLight));
    }

    PointDouble3DLight(const PointDouble3D& _other)
    {
        *this = _other;
    }

    PointDouble3DLight(const double& x, const double& y, const double& z)
    {
        this->x = x;
        this->y = y;
        this->z = static_cast<uint16_t>(z);
    }

    ~PointDouble3DLight() {}
    inline PointDouble3DLight& operator=(const PointDouble3DLight& point)
    {
        memcpy(this, &point, sizeof(PointDouble3DLight));
        return *this;
    }
    inline PointDouble3DLight& operator=(const PointDouble3D& point)
    {
        this->x = point.x;
        this->y = point.y;
        this->z = static_cast<uint16_t>(point.z);
        return *this;
    }
    inline bool operator==(const double& d) const
    {
        return (0. == std::abs(d - x) && 0. == std::abs(d - y) && 0. == std::abs(d - z));
    }
    inline bool operator!=(const double& d) const
    {
        return (0. != std::abs(d - x) || 0. != std::abs(d - y) || 0. != std::abs(d - z));
    }
    inline bool operator!=(const PointDouble3DLight& d)
    {
        return (0. != std::abs(d.x - x) || 0. != std::abs(d.y - y) || 0. != std::abs(d.z - z));
    }

    friend std::ostream &operator<<( std::ostream &output, const PointDouble3DLight &_p )
    {
      output << std::setprecision(6) << std::fixed << "X: " << _p.x << ", Y: " << _p.y << ", Z: " << _p.z;
      return output;
    }

    inline double abs3()
    {
        return sqrt(x*x+y*y+z*z);
    }
//    inline void normalize()
//    {
//        double absVal = abs3();
//        x /= absVal;
//        y /= absVal;
//        z /= absVal;
//    }

    double y;
    double x;
    uint16_t z;
}PACK_STRUCT;

bool operator==(const double& v, const PointDouble3DLight& pp);

bool operator!=(const double& v, const PointDouble3DLight& pp);

struct MachineId
{
    bool operator==(const MachineId& other) const
    {
        return std::strcmp(id, other.id) == 0;
    }
    char id[33] { 0 };
};

union ModelVersion
{
    ModelVersion()
    {
        version = 0;
    }

    ModelVersion(const uint64_t& _val)
    {
        version = _val;
    }

    ModelVersion(const uint32_t& _major, const uint32_t& _minor, const uint32_t& _patch)
    {
        major = _major;
        minor = _minor;
        patch = _patch;
    }

    ModelVersion& operator=(const ModelVersion& _other)
    {
        version = _other.version;
        return *this;
    }

    bool operator==(const ModelVersion& _other)
    {
        return (version == _other.version);
    }

    bool operator!=(const ModelVersion& _other)
    {
        return (version != _other.version);
    }

    bool is_version_match(const ModelVersion& _other) const
    {
        return (major == _other.major && minor == _other.minor);
    }

    std::string to_string() const
    {
        std::stringstream _ss;
        _ss << major << "." << minor << "." << patch;
        return _ss.str();
    }

    static ModelVersion from_string(const std::string& _str)
    {
        ModelVersion _model;
        auto _list_strs = split_str(_str, '.');
        if(_list_strs.size() >= 3)
        {
            _model.major = std::stoul(_list_strs[0]);
            _model.minor = std::stoul(_list_strs[1]);
            _model.patch = std::stoul(_list_strs[2]);
        }
        return _model;
    }

    friend std::ostream& operator<<(std::ostream &os, const ModelVersion& m)
    {
        os << m.to_string() ;
        return os;
    }

    uint64_t version;
    struct{
        uint64_t major:24;
        uint64_t minor:16;
        uint64_t patch:24;
    }PACK_STRUCT;
}PACK_STRUCT;


#endif
