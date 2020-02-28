cmake_minimum_required(VERSION 3.15)
project(EKG-HRV-Analyser)

set(CMAKE_CXX_STANDARD 17)

find_package(ALGLIB REQUIRED)

add_executable(EKG-HRV-Analyser main.cpp)

target_link_libraries(EKG-HRV-Analyser alglib)