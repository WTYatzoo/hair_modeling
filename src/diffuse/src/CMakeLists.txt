cmake_minimum_required(VERSION 2.8)

project(diffuse)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++0x")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

#make eigen available
set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3/")
include_directories( "${EIGEN3_INCLUDE_DIR}" )

#make spectra available
set(SPECTRA_INCLUDE_DIR "/usr/local/include/spectra-master/include/")
include_directories( "${SPECTRA_INCLUDE_DIR}" )

#make boost available
set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

find_package(Boost 1.58.0 COMPONENTS date_time filesystem system program_options REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  message(${Boost_LIBRARIES})
else()
  message("--[info] can not find Boost")
endif()

#make opencv available
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR};${CAMKE_MODULE_PATH}")
find_package(OpenCV REQUIRED)

#make eigen SuiteSparse available
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR};${CAMKE_MODULE_PATH}")
find_package(SuiteSparse)
include_directories(${SuiteSparse_INCLUDE_DIRS})

aux_source_directory(. src)
file(GLOB * header)

link_directories("/usr/lib/x86_64-linux-gnu/")

add_executable(diffuse ${src} ${header})
target_link_libraries(diffuse GLU glut GL ${Boost_LIBRARIES} ${SuiteSparse_LIBRARIES} ${OpenCV_LIBS})
