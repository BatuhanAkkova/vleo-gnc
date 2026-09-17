cmake_minimum_required(VERSION 3.16)
project(gnc_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ── Compiler Optimisation ──────────────────────────────────────────────
if(MSVC)
    add_compile_options($<$<NOT:$<CONFIG:Debug>>:/O2> /fp:fast /openmp)
    add_compile_definitions(_USE_MATH_DEFINES)
else()
    add_compile_options(-O3 -march=native -ffast-math -fopenmp)
    add_link_options(-fopenmp)
endif()

# ── Eigen (header-only) ───────────────────────────────────────────────
# Try system Eigen first, fall back to FetchContent (header-only)
find_package(Eigen3 3.4 QUIET)
if(NOT Eigen3_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0
        GIT_SHALLOW    TRUE
    )
    # Populate without building — Eigen's CMakeLists has broken math-lib tests on some MSVC
    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
    endif()
    # Create a header-only INTERFACE target
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${eigen_SOURCE_DIR}"
    )
endif()

# ── nlohmann/json (header-only, for JSON data loading) ────────────────
include(FetchContent)
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(json)



# ── Source files ──────────────────────────────────────────────────────
# Gather all .cpp under src/
file(GLOB_RECURSE GNC_SOURCES CONFIGURE_DEPENDS "src/*.cpp")

# Exclude main.cpp from library sources
list(FILTER GNC_SOURCES EXCLUDE REGEX "main\\.cpp$")

# Core library: STATIC if we have .cpp files, INTERFACE if header-only
if(GNC_SOURCES)
    add_library(gnc_core STATIC ${GNC_SOURCES})
    target_include_directories(gnc_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(gnc_core PUBLIC Eigen3::Eigen nlohmann_json::nlohmann_json)
else()
    add_library(gnc_core INTERFACE)
    target_include_directories(gnc_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(gnc_core INTERFACE Eigen3::Eigen nlohmann_json::nlohmann_json)
endif()



# Main executable
add_executable(gnc_sim src/main.cpp)
target_link_libraries(gnc_sim PRIVATE gnc_core)

# ── GoogleTest ────────────────────────────────────────────────────────
option(BUILD_TESTS "Build unit tests" ON)
if(BUILD_TESTS)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.14.0
        GIT_SHALLOW    TRUE
    )
    # Prevent GoogleTest from overriding compiler/linker options on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    enable_testing()

    file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "tests/*.cpp")
    add_executable(gnc_tests ${TEST_SOURCES})
    target_link_libraries(gnc_tests PRIVATE gnc_core GTest::gtest_main)
    target_include_directories(gnc_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

    include(GoogleTest)
    gtest_discover_tests(gnc_tests)
endif()

