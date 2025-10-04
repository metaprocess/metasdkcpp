#include "ElapsedTimer.h"
#include <sstream>
#include <iomanip>
#include <iostream>

ElapsedTimer::ElapsedTimer(TimeUnit time_unit)
    : m_time_unit(time_unit)
{
    m_start = std::chrono::high_resolution_clock::now();
}

std::string ElapsedTimer::restart(const std::string& _msg2)
{
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = now - m_start;
    double _time_double = duration.count(); // Time in seconds
    m_start = now; // Restart the timer

    auto _msg = _msg2;
    if (_msg.empty()) {
        std::stringstream _stream;
        _stream << "module" << std::setw(4) << m_counter_module;
        _msg = _stream.str();
    }
    std::stringstream _stream;
    std::string unit_str{" seconds."};

    switch (m_time_unit) {
        case TimeUnit::Milliseconds:
            _stream << std::setw(30) << _msg << " -->" << std::setw(7) << std::setprecision(3) << std::fixed << _time_double << unit_str;
            break;
        case TimeUnit::Microseconds:
            _stream << std::setw(30) << _msg << " -->" << std::setw(10) << std::setprecision(6) << std::fixed << _time_double << unit_str;
            break;
    }

    m_counter_module++;
    return _stream.str();
}

void ElapsedTimer::print_and_restart(const std::string& _msg)
{
    std::cerr << restart(_msg) << "\n";
}
