cmake_minimum_required(VERSION 2.8)
project(prjece650)

if (EXISTS "${CMAKE_SOURCE_DIR}/minisat")
  add_subdirectory(minisat)
else()
  message (FATAL_ERROR "MiniSat not found. Use `git clone https://git.uwaterloo.ca/ece650-1231/minisat.git` to check it out.")
endif()

# file(REMOVE ${BUILD_DIR}/example.csv)
# Copy example.csv to the build directory
# configure_file(${SOURCE_DIR}/example.csv ${BUILD_DIR}/example.csv COPYONLY)
# configure_file(example.csv ${CMAKE_CURRENT_BINARY_DIR}/example.csv COPYONLY)

## ECE650 Specific configuration start here
## Project specific configurations go here

# enable c++11 support
set (CMAKE_CXX_FLAGS "-std=c++11 -Wall ${CMAKE_CXX_FLAGS}")

# expose minisat header files
include_directories(${CMAKE_SOURCE_DIR}/minisat)

# create the main executable
## add additional .cpp files if needed
add_executable(prjece650 prjece650.cpp)
# link MiniSAT libraries
target_link_libraries(prjece650 minisat-lib-static)
target_link_libraries(prjece650 ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(prjece650 pthread)
target_link_libraries(prjece650 -lpthread)


# create the executable for tests
# add_executable(test test.cpp)
# link MiniSAT libraries
# target_link_libraries(test minisat-lib-static)
