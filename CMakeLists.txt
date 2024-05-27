cmake_minimum_required (VERSION 3.16)
project ("project")

# file (GLOB <variable> ... ): match files and set variable
file(GLOB SRCS "./src/*.cpp" "./src/*.c" "./imgui/*.cpp" "./imgui/backends/*.cpp")

add_executable(project ${SRCS})

target_include_directories(project PRIVATE "./include" "./imgui" "./imgui/backends")

set(LINUX FALSE)

if(LINUX)
    target_link_libraries(project dl GL glfw)
else()
    target_link_directories(project PRIVATE "./lib")
    target_link_libraries(project glfw3)
endif()
