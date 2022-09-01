#ifndef TIMER_HPP
#define TIMER_HPP

#include <iostream>
#include <vector>
#include <string>
#include <chrono>

class timer {
public:
    std::vector<double> times;
    std::vector<std::string> names;
    
    timer(std::vector<std::string> names) : times(names.size(), 0), names(names) { }
    
    int current_section = -1;
    std::chrono::time_point<std::chrono::high_resolution_clock> t_last;
    
    void start_section(int i) {
        auto t = std::chrono::high_resolution_clock::now();
        
        double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t - t_last).count();
        if (current_section >= 0)
            times[current_section] += dt;
        
        t_last = t;
        current_section = i;
    }
    
    void end_section() {
        auto t = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t - t_last).count();
        if (current_section >= 0)
            times[current_section] += dt;
        current_section = -1;
    }
    
    std::string format_times() {
        std::stringstream ss;
        ss << names[0] << ": " << times[0] * 1e-6 << "ms";
        for (int i = 1; i < (int)times.size(); i++)
            ss << ", " << names[i] << ": " << times[i] * 1e-6 << "ms";
        return ss.str();
    }
    
};

#endif


