#ifndef __ELAPSED_TIME_H__
#define __ELAPSED_TIME_H__

#include <chrono>
#include <string>

class ElapsedTimer {
public:
    enum class TimeUnit { Milliseconds, Microseconds };

    ElapsedTimer(TimeUnit time_unit = TimeUnit::Milliseconds);
    std::string restart(const std::string& _msg = "");
    void print_and_restart(const std::string& _msg = "");

private:
    std::chrono::high_resolution_clock::time_point m_start;
    int m_counter_module{1};
    TimeUnit m_time_unit;
};

#endif