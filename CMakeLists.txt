cmake_minimum_required(VERSION 3.16)
project(Seifenhaut LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets OpenGLWidgets)

if(MSVC)
    # M_PI and friends; keep windows.h from clobbering std::min/std::max.
    add_compile_definitions(_USE_MATH_DEFINES NOMINMAX)
endif()

add_executable(seifenhaut
    src/main.cpp
    src/MainWindow.cpp
    src/MainWindow.h
    src/FilmWidget.cpp
    src/FilmWidget.h
    src/SoapFilm.cpp
    src/SoapFilm.h
    src/Catenoid.h
)
target_link_libraries(seifenhaut PRIVATE Qt6::Widgets Qt6::OpenGLWidgets)

add_executable(verify tests/verify.cpp src/SoapFilm.cpp)
target_include_directories(verify PRIVATE src)
