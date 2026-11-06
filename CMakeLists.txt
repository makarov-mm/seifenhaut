cmake_minimum_required(VERSION 3.21)
project(Seifenhaut LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Qt6 REQUIRED COMPONENTS Widgets OpenGLWidgets)
qt_standard_project_setup()

qt_add_executable(seifenhaut
    src/main.cpp
    src/MainWindow.cpp
    src/MainWindow.h
    src/FilmWidget.cpp
    src/FilmWidget.h
    src/SoapFilm.cpp
    src/SoapFilm.h
    src/Catenoid.h
    src/Vec3.h
)
target_link_libraries(seifenhaut PRIVATE Qt6::Widgets Qt6::OpenGLWidgets)

if(MSVC)
    # the sources are UTF-8 (em dashes, Greek letters in the UI strings)
    target_compile_options(seifenhaut PRIVATE /utf-8 /W4)
else()
    target_compile_options(seifenhaut PRIVATE -Wall -Wextra)
endif()

set_target_properties(seifenhaut PROPERTIES
    WIN32_EXECUTABLE ON
    MACOSX_BUNDLE ON
)
