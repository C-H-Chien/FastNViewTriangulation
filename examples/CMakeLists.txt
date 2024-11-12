project(NViewsTrian-Examples CXX)


# Find Eigen library
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()


# Add the Eigen include directories
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(exp_gen ../utils/generatePointCloud.h ../utils/generatePointCloud.cpp)
#add_library(exp_ceres ../utils/ceresSolver.h ../utils/ceresSolver.cpp)

# Basic example
add_executable(example_base ${CMAKE_CURRENT_SOURCE_DIR}/example_base.cpp)
add_executable(mwv_triangulation_main ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

target_link_libraries(example_base 
                              NViewsTrian
                              exp_gen )         

target_link_libraries(mwv_triangulation_main 
                              NViewsTrian
                              exp_gen          )

