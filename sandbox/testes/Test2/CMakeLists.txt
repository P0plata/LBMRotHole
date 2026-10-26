cmake_minimum_required(VERSION 3.16)
project(LBM_Test_Local LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

include(FetchContent)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz)
FetchContent_MakeAvailable(json)

# Junta os ficheiros locais desta pasta de teste
add_executable(run_test
    main.cpp
    initializer.cpp
    custom_init.cpp
    geometry_builder.cpp
    vtk_exporter.cpp
    utils.cpp
    solver_d2q9.cpp
)

target_link_libraries(run_test PRIVATE OpenMP::OpenMP_CXX nlohmann_json::nlohmann_json)


# Se não especificares nada no comando cmake, assume Release (-O3 -march=native) por defeito
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

target_compile_options(run_test PRIVATE
    -Wall -Wextra
    $<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>
    $<$<CONFIG:Debug>:-O0 -g -fsanitize=address,undefined>
)
target_link_options(run_test PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address,undefined>
)