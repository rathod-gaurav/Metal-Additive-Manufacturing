cmake_minimum_required(VERSION 3.14)
project(calphad_client_example CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Download the two header-only dependencies once, into the build folder.
set(DEPS ${CMAKE_BINARY_DIR}/deps)
if(NOT EXISTS ${DEPS}/httplib.h)
  file(DOWNLOAD https://raw.githubusercontent.com/yhirose/cpp-httplib/v0.18.3/httplib.h
       ${DEPS}/httplib.h SHOW_PROGRESS)
endif()
if(NOT EXISTS ${DEPS}/json.hpp)
  file(DOWNLOAD https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp
       ${DEPS}/json.hpp SHOW_PROGRESS)
endif()

find_package(Threads REQUIRED)
add_executable(calphad_example example.cpp)
target_include_directories(calphad_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${DEPS})
target_link_libraries(calphad_example PRIVATE Threads::Threads)
if(WIN32)
  target_link_libraries(calphad_example PRIVATE ws2_32)
endif()

add_executable(kks_example kks_example.cpp)
target_include_directories(kks_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${DEPS})
target_link_libraries(kks_example PRIVATE Threads::Threads)
if(WIN32)
  target_link_libraries(kks_example PRIVATE ws2_32)
endif()
