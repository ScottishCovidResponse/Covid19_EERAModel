#pragma once

#include <chrono>

namespace EERAModel {

/**
 * @class SimpleTimer
 * @brief Simple timer for recording how long code runs take
 */
class SimpleTimer
{
public:
    /**
     * @brief Constructor
     */
    SimpleTimer() : startTime_(std::chrono::system_clock::now()) {}
    
    /**
     * @brief get elapsed time
     * 
     * @return Number of seconds since the timer was created
     */
    double elapsedTime() 
    { 
        auto stopTime = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime_);

        return elapsed.count() * 0.001;
    }

private:
    /**
     * @private
     * @brief Start time of the run
     */
    std::chrono::time_point<std::chrono::system_clock> startTime_;
};

} // namespace EERAModel
