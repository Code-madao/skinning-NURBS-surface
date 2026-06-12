cmake_minimum_required (VERSION 3.16)
project ("project")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release mode by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# file (GLOB <variable> ... ): match files and set variable
file(GLOB SRCS "./src/*.cpp" "./src/*.c" "./imgui/*.cpp" "./imgui/backends/*.cpp")

add_executable(project ${SRCS})

target_include_directories(project PRIVATE "./include" "E:/wheels/eigen-3.4.0" "./imgui" "./imgui/backends")

if(MSVC)
    target_compile_options(project PRIVATE /wd4477 /wd4819)
endif()

if(WIN32)
    target_link_directories(project PRIVATE "./lib")
    target_link_libraries(project glfw3)
elseif(UNIX)
    target_link_libraries(project dl GL glfw)
else()
    message(FATAL_ERROR "Unsupported platform")
endif()
