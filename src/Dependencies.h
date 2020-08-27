#pragma once

#include <string>

namespace EERAModel {

/**
 * @class DependencyVersions
 * @brief Class providing run-time access to project dependency versions
 * 
 * Implementation is populated by CMake's configure_file command
 */
class DependencyVersions {
public:
    /**
     * @brief Return the CMake version in use
     */
    static std::string CMakeVersion();

    /**
     * @brief Return the compiler id in use
     */
    static std::string CompilerId();

    /**
     * @brief Return the compiler version in use
     */
    static std::string CompilerVersion();

    /**
     * @brief Return the GSL version in use
     */
    static std::string GSLVersion();
};

} // namespace EERAModel
