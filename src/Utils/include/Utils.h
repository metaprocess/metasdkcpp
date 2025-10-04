#ifndef __METAUTILS_H__
#define __METAUTILS_H__

#include <complex>

#include <unistd.h>
#include "Definitions.h"
#include <sys/sysinfo.h>
#include <sys/types.h>
#include <iostream>
#include <zconf.h>
#include <netinet/in.h>
#include <functional>
#include <thread>

#define MISSIONS_ACTIONS_FILE "actions.csv"
#define MISSIONS_NO_SCENARIO "noScenario"

using type_callback_voidstar = std::function<void(void*)>;
using type_callback_string = std::function<void(std::string)>;
#define SCRIPTFILENAME "restartDefaultScript.sh"


class Utils
{
//    Q_OBJECT

public:
    static void set_app_name(const std::string& _name);
    static const std::string& get_app_name(void);
	static Utils* self;
	static int parseLine(char* line);

	static int getVal(const char* _S);
	static int parseLineWithoutKb(char* line);



    static int getValWithoutKb(const char* _S);
    static void setThreadName(const char* name);
	static long unsigned int systime,TotalIdle,totTime,utime;
	static int numProcessors;
	static double SysCpuPercent,ProcCpuPercent;
	/**
	 * @brief contain utility function for applicataion
	 * 
	 */
	Utils();
	
	static Utils* instance();

    static MachineId getMachineId();


    static std::string doubleToStdString(const double& val, const int& prec);
	

	static inline int64_t getTickCountUs()
    {
		struct timespec ts;
        clock_gettime ( CLOCK_REALTIME, &ts );
        return static_cast<int64_t>( ( ts.tv_sec * 1e6 ) + ( ts.tv_nsec / 1.e3 ) );
	}

	static inline int64_t getTickCountMs()
	{
		struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        return static_cast<int64_t>((ts.tv_sec * 1e3) + (ts.tv_nsec / 1.e6));
	}

   

	static void ecef2lla(double x, double y, double z, double& lat, double& lon, double& alt);
	static void lla2ecef(const double& lat, const double& lon, const double& alt, double& x, double& y, double& z);
	static void lla2ecefMatlab(const double& lat, const double& lon, const double& alt, double& x, double& y, double& z);
	static void lla2ecefMetaproc(const double& lat, const double& lon, const double& alt, double& x, double& y, double& z);

	static inline void miliSleep(int64_t t)
	{
        std::this_thread::sleep_for(std::chrono::milliseconds(t));
	}

    static inline void microSleep(int64_t t)
    {
        usleep(t);
    }
    
    template<class T>
    static inline T Pow2(const T& x)
	{
		return x * x;
	}

	template<class T>
	static inline T Pow3(const T& x)
	{
		return x * x * x;
	}

	template<class T>
	static inline T Pow4(const T& x)
	{
		return x * x * x * x;
	}

	static inline void Rotate ( double& x, double& y, double t )
	{
		double r0, t0;
		r0 = sqrt ( x * x + y * y );
		t0 = atan2 ( x, y );
		t0 += t ;
		x = r0 * sin ( t0 );
		y = r0 * cos ( t0 );
	}

	
    static void distance(double& r, double& theta_d, const double &lat1_d, const double &lon1_d, const double &lat2_d, const double &lon2_d);
    static void LL2RT(double& theta_r, double& r, const double* const lat1_r, const double* const lon1_r, const double* const lat2_r, const double* const lon2_r);
    static std::tuple<double, double> LL2RT(const double& lat1_r, const double& lon1_r, const double& lat2_r, const double& lon2_r);
    static void LL2RT(double& theta_r, double& r, const double& lat1_r, const double& lon1_r, const double& lat2_r, const double& lon2_r);
    static void LL2RT_DEG(double& theta_d, double& r, const double& lat1_d, const double& lon1_d, const double& lat2_d, const double &lon2_d);
    static void LL2RT_DEG(double& theta_d, double& r, const PointDouble& p1_d, const PointDouble& p2_d);
    static std::tuple<double, double> RT2LL(const double& r, const double& theta_r, const double& lat1_r, const double& lon1_r);
    static void RT2LL(const double& r, const double& theta_r, const double& lat1_r, const double& lon1_r, double& lat2_r, double& lon2_r);
    static void RT2LL_DEG(const double& r, const double& theta_d, const double& lat1_d, const double& lon1_d, double& lat2_d, double& lon2_d);
    static void RT2LL_DEG(const double& r, const double& theta_d, const PointDouble& p1_d, PointDouble &p2_d);
	static uint8_t calcXorCrc(void* data, int len)
	{
		uint8_t* cData = (uint8_t*)data;
		uint8_t crc = 0;
		for(int i = 0; i < len; i++)
		{
			crc ^= cData[i];
		}
		return crc;
    }

	static inline double getAzimuthDifference(double a, double b)
	{
		double diff = fabs(a-b);
		if(diff > M_PI)
			diff = 2*M_PI - diff;

		return diff;
	}
	static inline double getAzimuthDifferenceDeg(double a, double b)
	{
		double diff = fabs(a-b);
		if(diff > 180)
			diff = 360 - diff;

		return diff;
	}

    static inline double getAzimuthDifferenceDirectional(const double& az1, const double& az2)
        {

            double dd = az1 - az2;
            double diff = fabs(dd);
            if(diff > M_PI)
            {
                if(dd > 0)
                    dd -= 2*M_PI;
                else
                    dd += 2*M_PI;
            }
            return dd;
        }
    static inline double getAzimuthDifferenceDegDirectional(const double& az1, const double& az2)
    {

        double dd = az1 - az2;
        double diff = fabs(dd);
        if(diff > 180)
        {
            if(dd > 0)
                dd -= 2*180;
            else
                dd += 2*180;
        }
        return dd;
    }

	static double OffcenterAngleConvert ( double teta0 , double& r, double x0, double y0 )
	{
		double D =sqrt ( x0*x0+y0*y0 );
		double R=1;
		double gama=0, teta1 = 0, beta=0 ;
		double alpha=atan2 ( y0,x0 );
		
		if ( teta0 == alpha || ( teta0 == alpha+2*M_PI ) )
		{
			teta1 = alpha;
			r = R-D;
		}
		else if ( teta0 == ( alpha + M_PI ) || ( teta0 == ( alpha + 3*M_PI ) ) )
		{
			teta1 = alpha + M_PI;
			r = R+D;
		}
		else
		{
			double x = D*sin ( M_PI+alpha-teta0 );
			gama = asin ( x );
			beta = teta0 - alpha - gama;
			teta1 = alpha + beta;
			r=D*sin ( beta ) /sin ( gama );
		}
		if(teta1 >= 2*M_PI)
			teta1 -= 2*M_PI;
		else if(teta1 < 0)
			teta1 += 2*M_PI;
		return teta1;
	}
	template <class T>
	static inline void MinOfArray ( const T* array, int size, T& element, int& index )
		{
			if ( size < 1 )
			{
				index = -1;
				return;
			}
			else if ( size == 1 )
			{
				element = array[0];
				index = 0;
				return;
			}

			element = array[0];
			index = 0;
			for ( int i = 1; i < size; i++ )
			{
				if ( array[i] < element )
				{
					element = array[i];
					index = i;
				}
			}
		}
    template <class T>
    static inline void MaxOfArray ( const T* array, int size, T& element, int& index )
    {
        if ( size < 1 )
        {
            index = -1;
            return;
        }
        else if ( size == 1 )
        {
            element = array[0];
            index = 0;
            return;
        }

        element = array[0];
        index = 0;
        for ( int i = 1; i < size; i++ )
        {
            if ( array[i] > element )
            {
                element = array[i];
                index = i;
            }
        }
    }
    template <class T>
    static inline void MinMaxOfArray ( const T* array, int size, T& minElement, T& maxElement)
    {
        if ( size < 1 )
        {
            return;
        }
        else if ( size == 1 )
        {
            maxElement = minElement = array[0];
            return;
        }

        maxElement = minElement = array[0];
        for ( int i = 1; i < size; i++ )
        {
            if ( array[i] > maxElement )
            {
                maxElement = array[i];
            }
            if ( array[i] < minElement )
            {
                minElement = array[i];
            }
        }
    }
		

    static void generate_random_vectors(const size_t num_of_rand_vectors, std::vector<Byte> &vec) {
        vec.resize(num_of_rand_vectors);
        for (size_t j = 0; j < num_of_rand_vectors; ++j) {
            vec[j] = rand() % 255;
        }
    }
	static void hamming(int len, double* buffer);
	static void blackman(int len, double* buffer);
	static void kaiser(int len, double alpha, double* buffer);
	static int factorial(int n);
	static double mbesseli0(double x);
	static double kaiserWindow(double alpha, double x, int N);

	static void getJetColor ( const double& val2, double& red, double& green, double& blue )
	{
		double value = val2;
		if ( value > 1 )
			value = 1;

		double fourValue = 4 * value ;
		red   = std::min ( fourValue - 1.5, -fourValue + 4.5 );
		green = std::min ( fourValue - 0.5, -fourValue + 3.5 );
		blue  = std::min ( fourValue + 0.5, -fourValue + 2.5 );

		red   = std::max ( 0., red );
		green = std::max ( 0., green );
		blue  = std::max ( 0., blue );

		red   = std::min ( red, 1. );
		green = std::min ( green, 1. );
		blue  = std::min ( blue, 1. );
	}

    static void zlib_compress(const void *in_data, const size_t &in_data_size, std::vector<uint8_t> &out_buffer, const int &comp_level = -1);
    static int zlib_uncompress(const void *in_data, size_t in_data_size, std::vector<uint8_t> &out_buffer);
    static int zlib_uncompress_native(const void *in_data, size_t in_data_size, uint8_t *out_buffer, const size_t &max_out_size);

    static std::string path_join_multiple(const std::vector<std::string>& _lst);
    static void update_missions_structure_params(const char *_parent_address, const char* _mission_name, const char* _scenario_name, char* _mission_track_file_path, char* _mission_raw_log_path, char* _mission_imaging_tif_path, char *_mission_isar_png_path, char* _mission_actions_file_path);

    static const ModelVersion get_package_version();
    static void print_version();
    static const std::string timestamp_to_iso_datetime_string(const int64_t& _timestamp);
    static void startup();  // Entry point tasks
    static void cleanup();  // Exit point tasks
private:

	// NativeThread interface
protected:
	virtual void Execute(void *);
};

class MathUtils //: public QObject
{
	public:
        static void print_complex_array(const std::vector<std::complex<double> >&, const std::string& _name);
        static void print_array(const double* _buffer, const size_t& _size, const std::string& _name);
        static double asind(const double& ag);
        static double atand(const double& ag);
		static double sind(const double& ag);
        static double tand(const double& ag);
		static double cosd(const double& ag);
		static double  calcStdev(double* input, int size, double& mean, double& var);
		static double  calcStdevComplex(ComplexType* input, int size, double& mean, double& var);
		static void cross(const double a[3], const double b[3], double C[3]);
		static double norm(const double* vec, const int& len);
		static inline double db(const double& input)
		{
            return 20.*log10(input);
		}
		static double sign(const double& val)
		{
			return (val == 0 ? 0 : ( val > 0 ? 1 : -1 ) );
		}
		template <class TYPE_DATA, class TYPE_COEF>
        static inline void filter(const TYPE_COEF* b, const TYPE_DATA* _x_in, TYPE_DATA* _y_out, int _size_coef, int _size_x_in, double _denominator = 1.)
		{
			int Xdim=_size_x_in + _size_coef -1;
			TYPE_DATA* _vec_tmp = new TYPE_DATA[Xdim];
			bzero(_vec_tmp, sizeof(_vec_tmp[0]) * (_size_coef - 1));

			memcpy(&_vec_tmp[_size_coef - 1], _x_in, sizeof(_x_in[0]) * _size_x_in);

            if(1. != _denominator && 0. != _denominator)
            {
                for(int n = 0; n < _size_x_in ; n++)
                {
                    _y_out[n] = 0;
                    for(int j = 0; j < _size_coef; j++)
                    {
                        _y_out[n] += _vec_tmp[n + _size_coef - 1 - j] * b[j];
                    }
                    _y_out[n] /= _denominator;
                }

            }
            else if(0. != _denominator)
            {
                for(int n = 0; n < _size_x_in; n++)
                {
                    _y_out[n] = 0;
                    for(int j = 0; j < _size_coef; j++)
                    {
                        _y_out[n] += _vec_tmp[n + _size_coef - 1 - j] * b[j];
                    }
                }
            }
            else if(0. == _denominator)
            {
                std::cerr << "filter coefficient a is zero\n";
                abort();
            }

			delete [] _vec_tmp;
		}
		
		template <class T>
		static inline void fftShift(T* input, int len)
		{
			assertm( len > 2, " len > 2");
			T* temp = new T[len];
			if((len % 2) != 0)
			{
				memcpy(temp, &input[(len / 2) + 1 ], (len / 2) * sizeof(input[0]));
				memcpy(&temp[len / 2] , &input[0], sizeof(input[0]));
				memcpy(&temp[(len / 2) + 1], &input[1], (len / 2) * sizeof(input[0]));
				memcpy(input, temp, len * sizeof(input[0]));
			}
			else
			{
				memcpy(temp, &input[(len / 2)], (len / 2) * sizeof(input[0]) );
				memcpy(&temp[(len / 2)], input, (len / 2) *sizeof(input[0]));
				memcpy(input, temp, len * sizeof(input[0]));
			}
			delete [] temp;
        }

		template <class T>
		static inline void ifftShift(T* input, int len)
		{
			assertm( len > 2, " len > 2");
			T* temp = new T[len];
			if((len % 2) != 0)
			{
				memcpy(temp, &input[ len / 2 ], (len / 2)*sizeof(input[0]) );
				memcpy(&temp[len / 2] , &input[len-1], sizeof(input[0]) );
				memcpy(&temp[(len/2) + 1 ], &input[0], (len/2) *sizeof(input[0]));
				memcpy(input, temp, len * sizeof(input[0]));
			}
			else
			{
				memcpy(temp, &input[ (len/2) ], (len / 2)*sizeof(input[0]) );
				memcpy(&temp[(len/2) ], input, (len/2) *sizeof(input[0]));
				memcpy(input, temp, len * sizeof(input[0]));
			}
			delete [] temp;
		}
		
		static inline double meanArray(double* in, int len)
		{
			double mm = 0;
			for(int i = 0; i < len; i++)
			{
				mm += in[i];
			}
			mm /= (double)len;
			return mm;
		}
		static inline double absComplex(const double& re, const double& im)
		{
			return sqrt( re*re + im*im );
		}
		static inline double normalizeAngleRadian(const double& rad)
		{
			double dd = rad;
			return normalizeAngleRadian(dd);
		}
		static inline double normalizeAngleRadian(double& rad)
		{
			rad = fmod(rad, 2*M_PI);
			if(rad < 0)
				rad += 2*M_PI;
			return rad;
		}
		static inline double normalizeAngleDeg(const double& deg)
		{
			double dd = deg;
			return normalizeAngleDeg(dd);
		}
		static inline double normalizeAngleDeg(double& deg)
		{
			deg = fmod(deg, 360.);
			if(deg < 0)
				deg += 360.;
			return deg;
		}

		static void getJetColor ( const double& val, double& red, double& green, double& blue )
		{
			double value = val;
			if ( value > 1 )
				value = 1;

			double fourValue = 4 * value ;
			red   = std::min ( fourValue - 1.5, -fourValue + 4.5 )/* * 255.*/;
			green = std::min ( fourValue - 0.5, -fourValue + 3.5 )/* * 255.*/;
			blue  = std::min ( fourValue + 0.5, -fourValue + 2.5 )/* * 255.*/;

			red   = std::max ( 0., red );
			green = std::max ( 0., green );
			blue  = std::max ( 0., blue );

			red   = std::min ( red, 1. );
			green = std::min ( green, 1. );
			blue  = std::min ( blue, 1. );
		}



};


#endif // UTILS_H
