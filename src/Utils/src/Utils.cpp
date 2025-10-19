#include "Utils.h"
#include <iostream>
#include <assert.h>

#include <ifaddrs.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include "KaiserWindow.h"
#include <fstream>
#include <numeric>
#include <unistd.h>
#include <vector>
#include <zlib.h>
#include <sys/prctl.h>
#include "filesystem"
#include <sys/ioctl.h>
#include <net/if.h>
#include <arpa/nameser.h>
#include <resolv.h>
#include <bitset>
#include <regex>


namespace {
    // Helper class for automatic startup/cleanup
    struct StartupCleanupHandler {
        StartupCleanupHandler() { Utils::startup(); }
        ~StartupCleanupHandler() { Utils::cleanup(); }
    };

    // Static instance ensures constructor/destructor are called
    StartupCleanupHandler globalHandler;
}

long unsigned int Utils::systime=0;
long unsigned int  Utils::TotalIdle=0;
long unsigned int  Utils::totTime=0;
long unsigned int  Utils::utime=0;
Utils* Utils::self = nullptr;
int Utils::numProcessors=0;

double Utils::SysCpuPercent=0;
double Utils::ProcCpuPercent=0;

std::string g_str_app_name;

std::vector<size_t> get_cpu_times() {
    std::ifstream proc_stat("/proc/stat");
    proc_stat.ignore(5, ' '); // Skip the 'cpu' prefix.
    std::vector<size_t> times;
    for (size_t time; proc_stat >> time; times.push_back(time));
    proc_stat.close();
    return times;
}

bool get_cpu_times(size_t &idle_time, size_t &total_time) {
    const std::vector<size_t> cpu_times = get_cpu_times();
    if (cpu_times.size() < 4)
        return false;
    idle_time = cpu_times[3];
    total_time = std::accumulate(cpu_times.begin(), cpu_times.end(), 0);
    return true;
}

void cpuUtilizationThread()
{
    size_t previous_idle_time=0, previous_total_time=0;
    for (size_t idle_time, total_time; get_cpu_times(idle_time, total_time); sleep(1)) {
        const float idle_time_delta = idle_time - previous_idle_time;
        const float total_time_delta = total_time - previous_total_time;
        const float utilization = 100.0 * (1.0 - idle_time_delta / total_time_delta);
        Utils::SysCpuPercent = utilization;
        //            std::cout << utilization << '%' << std::endl;
        previous_idle_time = idle_time;
        previous_total_time = total_time;
    }
}

void Utils::LL2RT(double& theta_r, double& r, const double * const lat1_r, const double * const lon1_r, const double * const lat2_r, const double * const lon2_r)
{
    double  A = 0, B = 0, C = 0 , L = 0, U1 = 0 , U2 = 0 , c2_sgm = 0 , c_sg = 0  ;
    double c_sqa = 0, d_sg = 0, itr = 0 , landa = 0 , landap = 0 , reverse_Azmt = 0 , sgm = 0, s_alp = 0 , s_sgm = 0 , uSq = 0  ;

    L = ( *lon2_r - *lon1_r );
    U1 = atan ( ( 1 - wgs84Meter.Flattening ) * tan ( *lat1_r ) );
    U2 = atan ( ( 1 - wgs84Meter.Flattening ) * tan ( *lat2_r ) );

    landa = L;
    itr = 100;
    landap = 100;
    while ( fabs ( landa - landap ) > 1e-12 && itr > 0 )
    {
        s_sgm = sqrt ( ( cos ( U2 ) * sin ( landa ) ) * ( cos ( U2 ) * sin ( landa ) ) + ( cos ( U1 ) * sin ( U2 ) - sin ( U1 ) * cos ( U2 ) * cos ( landa ) ) * ( cos ( U1 ) * sin ( U2 ) - sin ( U1 ) * cos ( U2 ) * cos ( landa ) ) );
        c_sg = sin ( U1 ) * sin ( U2 ) + cos ( U1 ) * cos ( U2 ) * cos ( landa );
        sgm = atan2 ( s_sgm, c_sg );
        s_alp = cos ( U1 ) * cos ( U2 ) * sin ( landa ) / s_sgm;
        c_sqa = 1 - s_alp * s_alp;
        c2_sgm = c_sg - 2 * sin ( U1 ) * sin ( U2 ) / c_sqa;
        C = wgs84Meter.Flattening / 16 * c_sqa * ( 4 + wgs84Meter.Flattening * ( 4 - 3 * c_sqa ) );
        landap = landa;
        landa = L + ( 1 - C ) * wgs84Meter.Flattening * s_alp * ( sgm + C * s_sgm * ( c2_sgm + C * c_sg * ( -1 + 2 * c2_sgm * c2_sgm ) ) );
        itr = itr - 1;
    }

    uSq = c_sqa * wgs84Meter.SquaredEccentricity;
    A = 1 + uSq / 16384 * ( 4096 + uSq * ( -768 + uSq * ( 320 - 175 * uSq ) ) );
    B = uSq / 1024 * ( 256 + uSq * ( -128 + uSq * ( 74 - 47 * uSq ) ) );
    d_sg = B * s_sgm * ( c2_sgm + B / 4 * ( c_sg * ( -1 + 2 * c2_sgm * c2_sgm ) - B / 6 * c2_sgm * ( -3 + 4 * s_sgm * s_sgm ) * ( -3 + 4 * c2_sgm * c2_sgm ) ) );
    r = wgs84Meter.SemiminorAxis * A * ( sgm - d_sg );
    //R = R.toFixed(3);% // round to 1mm precision

    //note: to return initial/final bearings in addition to distance, use something like:
    theta_r = atan2 ( cos ( U2 ) * sin ( landa ),  cos ( U1 ) * sin ( U2 ) - sin ( U1 ) * cos ( U2 ) * cos ( landa ) );
    reverse_Azmt = atan2 ( cos ( U1 ) * sin ( landa ), -sin ( U1 ) * cos ( U2 ) + cos ( U1 ) * sin ( U2 ) * cos ( landa ) );
}

std::tuple<double, double> Utils::LL2RT(const double &lat1_r, const double &lon1_r, const double &lat2_r, const double &lon2_r)
{
    double lt1 = lat1_r;
    double ln1 = lon1_r;
    double lt2 = lat2_r;
    double ln2 = lon2_r;

    double theta_r{0}, r{0};
    LL2RT(theta_r,r, &lt1, &ln1, &lt2, &ln2);
    return {r, theta_r};
}

void Utils::distance(double &r, double &theta_d, const double& lat1_d, const double& lon1_d, const double& lat2_d, const double& lon2_d)
{
    double lat1_r = lat1_d * DEGREE_TO_RADIAN_UNIT;
    double lon1_r = lon1_d * DEGREE_TO_RADIAN_UNIT;
    double lat2_r = lat2_d * DEGREE_TO_RADIAN_UNIT;
    double lon2_r = lon2_d * DEGREE_TO_RADIAN_UNIT;
    double  A = 0, B = 0, C = 0 , L = 0, U1 = 0 , U2 = 0 , c2_sgm = 0 , c_sg = 0  ;
    double c_sqa = 0, d_sg = 0, itr = 0 , landa = 0 , landap = 0 , reverse_Azmt = 0 , sgm = 0, s_alp = 0 , s_sgm = 0 , uSq = 0  ;

    L = ( lon2_r - lon1_r );
    U1 = atan ( ( 1 - wgs84Meter.Flattening ) * tan ( lat1_r ) );
    U2 = atan ( ( 1 - wgs84Meter.Flattening ) * tan ( lat2_r ) );

    landa = L;
    itr = 100;
    landap = 100;
    while ( fabs ( landa - landap ) > 1e-12 && itr > 0 )
    {
        s_sgm = sqrt ( ( cos ( U2 ) * sin ( landa ) ) * ( cos ( U2 ) * sin ( landa ) ) + ( cos ( U1 ) * sin ( U2 ) - sin ( U1 ) * cos ( U2 ) * cos ( landa ) ) * ( cos ( U1 ) * sin ( U2 ) - sin ( U1 ) * cos ( U2 ) * cos ( landa ) ) );
        c_sg = sin ( U1 ) * sin ( U2 ) + cos ( U1 ) * cos ( U2 ) * cos ( landa );
        sgm = atan2 ( s_sgm, c_sg );
        s_alp = cos ( U1 ) * cos ( U2 ) * sin ( landa ) / s_sgm;
        c_sqa = 1 - s_alp * s_alp;
        c2_sgm = c_sg - 2 * sin ( U1 ) * sin ( U2 ) / c_sqa;
        C = wgs84Meter.Flattening / 16 * c_sqa * ( 4 + wgs84Meter.Flattening * ( 4 - 3 * c_sqa ) );
        landap = landa;
        landa = L + ( 1 - C ) * wgs84Meter.Flattening * s_alp * ( sgm + C * s_sgm * ( c2_sgm + C * c_sg * ( -1 + 2 * c2_sgm * c2_sgm ) ) );
        itr = itr - 1;
    }

    uSq = c_sqa * wgs84Meter.SquaredEccentricity;
    A = 1 + uSq / 16384 * ( 4096 + uSq * ( -768 + uSq * ( 320 - 175 * uSq ) ) );
    B = uSq / 1024 * ( 256 + uSq * ( -128 + uSq * ( 74 - 47 * uSq ) ) );
    d_sg = B * s_sgm * ( c2_sgm + B / 4 * ( c_sg * ( -1 + 2 * c2_sgm * c2_sgm ) - B / 6 * c2_sgm * ( -3 + 4 * s_sgm * s_sgm ) * ( -3 + 4 * c2_sgm * c2_sgm ) ) );
    r = wgs84Meter.SemiminorAxis * A * ( sgm - d_sg );
    //R = R.toFixed(3);% // round to 1mm precision

    //note: to return initial/final bearings in addition to distance, use something like:
    theta_d = atan2 ( cos ( U2 ) * sin ( landa ),  cos ( U1 ) * sin ( U2 ) - sin ( U1 ) * cos ( U2 ) * cos ( landa ) );
    reverse_Azmt = atan2 ( cos ( U1 ) * sin ( landa ), -sin ( U1 ) * cos ( U2 ) + cos ( U1 ) * sin ( U2 ) * cos ( landa ) );

    theta_d *= RADIAN_TO_DEGREE_UNIT;

    if (theta_d < 0.0)
    {
        theta_d += 360.0;
    }
}

void Utils::LL2RT(double &theta_r, double &r, const double &lat1_r, const double &lon1_r, const double &lat2_r, const double &lon2_r)
{
    double lt1 = lat1_r;
    double ln1 = lon1_r;
    double lt2 = lat2_r;
    double ln2 = lon2_r;

    LL2RT(theta_r,r, &lt1, &ln1, &lt2, &ln2);
}

void Utils::LL2RT_DEG(double &theta_d, double &r, const double &lat1_d, const double &lon1_d, const double &lat2_d, const double& lon2_d)
{
    double _theta_r;
    LL2RT(_theta_r, r, lat1_d * DEGREE_TO_RADIAN_UNIT, lon1_d * DEGREE_TO_RADIAN_UNIT, lat2_d * DEGREE_TO_RADIAN_UNIT, lon2_d * DEGREE_TO_RADIAN_UNIT);
    theta_d = _theta_r * RADIAN_TO_DEGREE_UNIT;
}

void Utils::LL2RT_DEG(double &theta_d, double &r, const PointDouble &p1_d, const PointDouble &p2_d)
{
    double _theta_r;
    LL2RT(_theta_r, r, p1_d.y * DEGREE_TO_RADIAN_UNIT, p1_d.x * DEGREE_TO_RADIAN_UNIT, p2_d.y * DEGREE_TO_RADIAN_UNIT, p2_d.x * DEGREE_TO_RADIAN_UNIT);
    theta_d = MathUtils::normalizeAngleDeg(_theta_r * RADIAN_TO_DEGREE_UNIT);
}

std::tuple<double, double> Utils::RT2LL(const double &r, const double &theta_r, const double &lat1_r, const double &lon1_r)
{
    double _lat{0}, _lon{0};
    RT2LL(r, theta_r, lat1_r, lon1_r, _lat, _lon);
    return {_lat, _lon};
}

void Utils::RT2LL(const double& r, const double& theta_r, const double& lat1_r, const double& lon1_r, double& lat2_r, double& lon2_r)
{
    double aa = 0, A = 0, B = 0, C = 0, L = 0, c2sg = 0, ca1 = 0, /*csg0 = 0,*/ cs45 = 0, cu1 = 0, dsg = 0, lda = 0;
    double /* fbear = 0, */ s = 0, sgm_ = 0, sgm_1 = 0, s__p = 0, s_ap = 0 , s_ap1 = 0, s_sg = 0, s_U = 0, t_U1 = 0, t_M_P = 0, u0 = 0, csg = 0;

    s = r;
    aa = theta_r;
    s_ap1 = sin ( aa );
    ca1 = cos ( aa );

    t_U1 = ( 1.0 - wgs84Meter.Flattening ) * tan ( lat1_r );
    cu1 = 1.0 / sqrt ( ( 1.0 + t_U1 * t_U1 ) );
    s_U = t_U1 * cu1;
    sgm_1 = atan2 ( t_U1, ca1 );
    s_ap = cu1 * s_ap1;
    cs45 = 1.0 - s_ap * s_ap;
    u0 = cs45 * wgs84Meter.SquaredEccentricity;
    A = 1.0 + u0 / 16384.0 * ( 4096.0 + u0 * ( -768.0 + u0 * ( 320.0 - 175.0 * u0 ) ) );
    B = u0 / 1024.0 * ( 256.0 + u0 * ( -128.0 + u0 * ( 74.0 - 47.0 * u0 ) ) );

    sgm_ = s / ( wgs84Meter.SemiminorAxis * A );
    s__p = 2.0 * M_PI;
    int iter = 100;
    while ( fabs ( sgm_ - s__p ) > 1e-12  && iter > 0)
    {
        c2sg = cos ( 2.0 * sgm_1 + sgm_ );
        s_sg = sin ( sgm_ );
        csg = cos ( sgm_ );
        dsg = B * s_sg * ( c2sg + B / 4.0 * ( csg * ( -1.0 + 2.0 * c2sg * c2sg ) -  B / 6.0 * c2sg * ( -3.0 + 4.0 * s_sg * s_sg ) * ( -3.0 + 4.0 * c2sg * c2sg ) ) );
        s__p = sgm_;
        sgm_ = s / ( wgs84Meter.SemiminorAxis * A ) + dsg;
        iter--;
    }

    t_M_P = s_U * s_sg - cu1 * csg * ca1;
    lat2_r = atan2 ( s_U * csg + cu1 * s_sg * ca1, ( 1.0 - wgs84Meter.Flattening ) * sqrt ( s_ap * s_ap + t_M_P * t_M_P ) );
    lda = atan2 ( s_sg * s_ap1, cu1 * csg - s_U * s_sg * ca1 );
    C = wgs84Meter.Flattening / 16.0 * cs45 * ( 4.0 + wgs84Meter.Flattening * ( 4.0 - 3.0 * cs45 ) );
    L = lda - ( 1.0 - C ) * wgs84Meter.Flattening * s_ap * ( sgm_ + C * s_sg * ( c2sg + C * csg * ( -1 + 2 * c2sg * c2sg ) ) );

    //hh=fmod( ((*lon1)+L+3*M_PI),(2*M_PI) ) - M_PI;
    lon2_r = fmod ( ( ( lon1_r ) + L + 3.0 * M_PI ) , 2.0 * M_PI ) - M_PI; // normalise to -180...+180
    //
    // fbear = atan2 ( s_ap, -t_M_P ); // final bearing, if required
}

void Utils::RT2LL_DEG(const double &r, const double &theta_d, const double &lat1_d, const double &lon1_d, double& lat2_d, double& lon2_d)
{
    double lat_rad = 0, lon_rad = 0;
    RT2LL(r, theta_d * DEGREE_TO_RADIAN_UNIT, lat1_d * DEGREE_TO_RADIAN_UNIT, lon1_d * DEGREE_TO_RADIAN_UNIT, lat_rad, lon_rad);
    lat2_d = lat_rad * RADIAN_TO_DEGREE_UNIT;
    lon2_d = lon_rad * RADIAN_TO_DEGREE_UNIT;
}

void Utils::RT2LL_DEG(const double &r, const double &theta_d, const PointDouble &p1_d, PointDouble &p2_d)
{
    double lat_rad = 0, lon_rad = 0;
    RT2LL(r, theta_d * DEGREE_TO_RADIAN_UNIT, p1_d.y * DEGREE_TO_RADIAN_UNIT, p1_d.x * DEGREE_TO_RADIAN_UNIT, lat_rad, lon_rad);
    p2_d.y = lat_rad * RADIAN_TO_DEGREE_UNIT;
    p2_d.x = lon_rad * RADIAN_TO_DEGREE_UNIT;
}

Utils::Utils()
{
}

Utils *Utils::instance()
{
    if(nullptr == self)
    {
        self = new Utils();
    }

    return self;
}

MachineId Utils::getMachineId()
{
    std::ifstream file("/etc/machine-id");
    std::string machineID;

    if (file.is_open()) {
        std::getline(file, machineID);
        file.close();
    }


    MachineId id{ };
    memcpy(id.id, machineID.c_str(), sizeof(id.id));
    id.id[32] = '\0';
    return id;
}

std::string Utils::doubleToStdString(const double& val, const int& prec)
{
    std::ostringstream ss1;
    ss1.precision(prec);
    ss1 << std::fixed << val;
    std::string _str = ss1.str();
    return _str;
}
int Utils::getVal(const char* _S)
{
    //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    int len= strlen(_S);
    if(file)
    {
        while (fgets(line, 128, file) != nullptr)
        {
            if (strncmp(line, _S, len) == 0)
            {
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
    }
    return result;
}

int Utils::getValWithoutKb(const char* _S)
{
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    if(file)
    {
        int len= strlen(_S);

        while (fgets(line, 128, file) != nullptr)
        {
            if (strncmp(line, _S, len) == 0)
            {
                result = parseLineWithoutKb(line);
                break;
            }
        }
        fclose(file);
    }
    return result;
}

void Utils::setThreadName(const char *name)
{
    prctl(PR_SET_NAME,name,0,0,0);
}


int Utils::parseLine(char* line)
{
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}
int Utils::parseLineWithoutKb(char* line)
{
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    i = atoi(line);
    return i;
}





void MathUtils::print_complex_array(const std::vector<std::complex<double>>& _list_complex, const std::string& _name)
{
    // Set the output format to fixed and the precision to 8
    // std::cout << std::fixed << std::setprecision(8);

    int _ctr = 0;
    std::cout << "<" << _name << ">=\n";
    for(const auto& c : _list_complex) {
        std::cout <<
         std::setfill('0') << 
         std::setprecision(0) << std::setw(4) <<_ctr++ << 
         std::setw(0) << ".(" <<
         std::setw(15) << std::setprecision(8) << c.real() << 
         std::setw(0) << ", " << 
         std::setw(15) << c.imag() << 
         std::setw(0) << ") ";
        if(0 == (_ctr % 4))
        {
            std::cerr << "\n";
        }
    }

    std::cout << std::endl<<"</" << _name << ">\n"; // New line after printing the array
}

void MathUtils::print_array(const double *_buffer, const size_t &_size, const std::string &_name)
{
    // Set the output format to fixed and the precision to 8
    // std::cout << std::fixed << std::setprecision(8);

    std::cout << "<" << _name << ">=\n";
    for(size_t _idx = 0; _idx < _size; _idx++) {
        std::cout <<
         std::setfill(' ') << 
         std::setprecision(0) << std::setw(4) <<_idx << 
         std::setw(0) << ":(" <<
         std::setw(15) << std::setprecision(8) << _buffer[_idx] <<  
         std::setw(0) << ") ";
        if(0 == (_idx % 4))
        {
            std::cerr << "\n";
        }
    }

    std::cout << std::endl<<"</" << _name << ">\n"; // New line after printing the array
}

double MathUtils::asind(const double& ag)
{
    return asin(ag) * RADIAN_TO_DEGREE_UNIT;
}

double MathUtils::atand(const double& ag)
{
    return atan(ag) * RADIAN_TO_DEGREE_UNIT;
}

double MathUtils::sind(const double &ag)
{
    return sin(ag * DEGREE_TO_RADIAN_UNIT);
}

double MathUtils::tand(const double &ag)
{
    return tan(ag * DEGREE_TO_RADIAN_UNIT);
}

double MathUtils::cosd(const double &ag)
{
    return cos(ag * DEGREE_TO_RADIAN_UNIT);
}

double  MathUtils::calcStdev(double* input, int size, double& mean, double& var)
{
    assert(size > 1);
    var = 0;
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
        sum += input[i];
    }

    mean = sum / (double)size;
    double sum2 = 0;
    for(int i = 0; i < size; i++)
    {
        sum2 += Utils::Pow2<double>(input[i] - mean);
    }
    var = sum2 / (double)(size-1);
    double stdev = sqrt( var );
    return stdev;
}

double MathUtils::calcStdevComplex(ComplexType *input, int size, double &mean, double &var)
{
    assert(size > 1);
    var = 0;
    double sum = 0;
    for(int i = 0; i < size; i++)
    {
        sum += input[i].abs();
    }

    mean = sum / (double)size;
    double sum2 = 0;
    for(int i = 0; i < size; i++)
    {
        sum2 += Utils::Pow2<double>(input[i].abs() - mean);
    }
    var = sum2 / (double)(size-1);
    double stdev = sqrt( var );
    return stdev;

}

void MathUtils::cross(const double a[], const double b[], double C[])
{
    C[0] = a[1]*b[2]-a[2]*b[1];
    C[1] = a[2]*b[0]-a[0]*b[2];
    C[2] = a[0]*b[1]-a[1]*b[0];
}

double MathUtils::norm(const double *vec, const int &len)
{
    double sum = 0;
    for(int i = 0; i < len; i++)
    {
        sum += Utils::Pow2(vec[i]);
    }
    sum /= (double)len;
    return sum;
}

#define CHUNK (4*1024)

void Utils::zlib_compress(const void *in_data, const size_t& in_data_size, std::vector<uint8_t> &out_buffer, const int& comp_level)
{
//    std::vector<uint8_t> buffer;

    const size_t BUFSIZE = CHUNK;
    uint8_t temp_buffer[BUFSIZE];

    z_stream strm;
    strm.zalloc = 0;
    strm.zfree = 0;
    strm.next_in = (Bytef*)(in_data);
    strm.avail_in = in_data_size;
    strm.next_out = temp_buffer;
    strm.avail_out = BUFSIZE;

    deflateInit(&strm, comp_level);

    while (strm.avail_in != 0)
    {
        int res = deflate(&strm, Z_NO_FLUSH);
        assert(res == Z_OK);
        if (strm.avail_out == 0)
        {
            out_buffer.insert(out_buffer.end(), temp_buffer, temp_buffer + BUFSIZE);
            strm.next_out = temp_buffer;
            strm.avail_out = BUFSIZE;
        }
    }

    int deflate_res = Z_OK;
    while (deflate_res == Z_OK)
    {
        if (strm.avail_out == 0)
        {
            out_buffer.insert(out_buffer.end(), temp_buffer, temp_buffer + BUFSIZE);
            strm.next_out = temp_buffer;
            strm.avail_out = BUFSIZE;
        }
        deflate_res = deflate(&strm, Z_FINISH);
    }

    assert(deflate_res == Z_STREAM_END);
    out_buffer.insert(out_buffer.end(), temp_buffer, temp_buffer + BUFSIZE - strm.avail_out);
    deflateEnd(&strm);

    //    out_data.swap(out_buffer);
}

int Utils::zlib_uncompress(const void *in_data, size_t in_data_size2, std::vector<uint8_t> &out_buffer)
{

    size_t in_data_size = in_data_size2;
    int ret = 0;
    unsigned have = 0;
    z_stream strm;
//    unsigned char in[CHUNK];
    int64_t in_data_index = 0;
    unsigned char out[CHUNK];
    /* allocate inflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK)
        return ret;
    /* decompress until deflate stream ends or end of file */
    do {
        size_t min_in_data_size = std::min(static_cast<size_t>(CHUNK), in_data_size);
        strm.avail_in = min_in_data_size;
//        if (ferror(source)) {
//            (void)inflateEnd(&strm);
//            return Z_ERRNO;
//        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = (Bytef*)(in_data) + in_data_index;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
                case Z_NEED_DICT:
                    ret = Z_DATA_ERROR;     /* and fall through */
                case Z_DATA_ERROR:
                case Z_MEM_ERROR:
                    (void)inflateEnd(&strm);
                    return ret;
            }
            have = CHUNK - strm.avail_out;
            out_buffer.insert(out_buffer.end(), out, out+have);

//            if (fwrite(out, 1, have, dest) != have || ferror(dest)) {
//                (void)inflateEnd(&strm);
//                return Z_ERRNO;
//            }

        } while (strm.avail_out == 0);
        in_data_index += min_in_data_size;
        in_data_size -= min_in_data_size;

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    /* clean up and return */
    (void)inflateEnd(&strm);
    return (ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR);

}

int Utils::zlib_uncompress_native(const void *in_data, size_t in_data_size2, uint8_t* out_buffer, const size_t& max_out_size)
{
    size_t in_data_size = in_data_size2;
    assert(out_buffer != nullptr);
    int ret = 0;
    unsigned have = 0;
    z_stream strm;
//    unsigned char in[CHUNK];
    int64_t in_data_index = 0;
    uint8_t out[CHUNK];
    /* allocate inflate state */
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    size_t out_buffer_size = 0;
    ret = inflateInit(&strm);
    if (ret != Z_OK)
        return ret;
    /* decompress until deflate stream ends or end of file */
    do {
        size_t min_in_data_size = std::min(static_cast<size_t>(CHUNK), in_data_size);
        strm.avail_in = min_in_data_size;
//        if (ferror(source)) {
//            (void)inflateEnd(&strm);
//            return Z_ERRNO;
//        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = (Bytef*)(in_data) + in_data_index;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
                case Z_NEED_DICT:
                    ret = Z_DATA_ERROR;     /* and fall through */
                case Z_DATA_ERROR:
                case Z_MEM_ERROR:
                    (void)inflateEnd(&strm);
                    return ret;
            }
            have = CHUNK - strm.avail_out;
//            out_buffer = reinterpret_cast<uint8_t*>(realloc(out_buffer, out_buffer_size + have));
            if((out_buffer_size + have) > max_out_size )
            {
                out_buffer_size = 0;
            }
            else
            {
                memcpy(out_buffer + out_buffer_size, out, have);
                out_buffer_size += have;
            }
//            out_buffer.add_to_end(out, have);
//            out_size = out_buffer.size();
//            out_buffer.insert(out_buffer.end(), out, out+have);

//            if (fwrite(out, 1, have, dest) != have || ferror(dest)) {
//                (void)inflateEnd(&strm);
//                return Z_ERRNO;
//            }

        } while (strm.avail_out == 0);
        in_data_index += min_in_data_size;
        in_data_size -= min_in_data_size;

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    /* clean up and return */
    (void)inflateEnd(&strm);
    return (ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR);

}

std::string Utils::path_join_multiple(const std::vector<std::string> &_lst)
{
    std::filesystem::path _path;
    int _ctr = 0;
    for(const std::string& _str: _lst)
    {
        auto _str_copy = _str;
        if(_ctr > 0 && std::filesystem::path::preferred_separator == _str[0])
        {
            _str_copy.erase(_str_copy.begin());
        }
        if(std::filesystem::path::preferred_separator == *(_str.end()-1))
        {
            _str_copy.pop_back();
        }
        _path.append(_str_copy);
        _ctr++;
    }
    return _path.string();
}

const ModelVersion Utils::get_package_version()
{
    return ModelVersion::from_string(PRJ_PACKAGE_VER);
}

void Utils::print_version()
{
    auto _obj = get_package_version();
    std::cerr << _obj << "\n";
}

const std::string Utils::timestamp_to_iso_datetime_string(const int64_t &_timestamp)
{
    time_t now = _timestamp;
//    time(&now);
    char buf[sizeof "2011-10-08T07:07:09Z"];
    strftime(buf, sizeof(buf), "%FT%TZ", gmtime(&now));
    return std::string(buf);
}

StdThreadPool &Utils::get_std_threadpool()
{
    static StdThreadPool _std_thread_pool;
    return _std_thread_pool;
}

void Utils::set_app_name(const std::string &_name)
{
    g_str_app_name = _name;
}

const std::string &Utils::get_app_name(void)
{
    return g_str_app_name;
}

void Utils::hamming(int len, double *buffer)
{
    for(int i = 0; i < len; i++)
    {
        buffer[i] = 0.54-.46*cos(2*M_PI*i/(double)(len-1) );
    }
}

void Utils::blackman(int len, double *buffer)
{
    for(int i = 0; i < len; i++)
    {
        buffer[i] = 0.42-.5*cos(2.*M_PI*(double)i/(double)(len-1) ) + 0.08 * cos( (4. * M_PI * (double)i) / double(len - 1) );
    }
}

void Utils::kaiser(int len, double alpha, double *buffer)
{
    KaiserWindow::buildWindow(len, alpha, buffer);
    //    for(int i = 0; i < len; i++)
    //    {
    //        buffer[i] = kaiserWindow(alpha, (double)i, len);
    //    }
}

int Utils::factorial(int n)
{
    if( n == 0)
    {
        return 1;
    }
    int val = 1;
    for( int idx = 1; idx <=n; ++idx)
    {
        val *= idx;
    }
    return val;
}

double Utils::mbesseli0(double x)
{
    const int MaxIter = 200;

    double sum = 0.0;
    for(int m = 0; m < MaxIter; ++m)
    {
        int factM = factorial(m); // m!
        double inc = Pow2<double>(1.0/factM * pow(x * 0.5, m)); // ( 1/(m!) * (x/2)^m )^2
        double frac = inc / sum;
        sum += inc;
        if( frac < 0.001)
        {
            break;
        }
    }
    return sum;

}

double Utils::kaiserWindow(double alpha, double x, int N)
{
    double b = M_PI * alpha * sqrt( 1.0 - pow( 2.0*x/(N-1.0) - 1.0, 2.0 ));
    double w = mbesseli0(b) / mbesseli0( M_PI * alpha);
    return w;
}

void Utils::Execute(void *)
{
    cpuUtilizationThread();
}


void Utils::startup()
{
    std::cerr << "GLOBAL STARTUP ROUTINE...\n";
}

void Utils::cleanup()
{
    std::cerr << "... GLOBAL CLEANUP ROUTINE\n";
}
