cmake_minimum_required(VERSION 2.8)

project(data_aug)

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

#make assimp availabl
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR};${CAMKE_MODULE_PATH}")
find_package(ASSIMP REQUIRED)
if(ASSIMP_FOUND)
    include_directories(${ASSIMP_INCLUDE_DIR})
endif() 

#make opencv available
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR};${CAMKE_MODULE_PATH}")
find_package(OpenCV REQUIRED)

#make openmesh available
if( (NOT DEFINED WALBERLA_BUILD_WITH_OPENMESH) OR WALBERLA_BUILD_WITH_OPENMESH )
   find_package( OpenMesh )
   if( OPENMESH_FOUND )
      set( WALBERLA_BUILD_WITH_OPENMESH ON CACHE BOOL "Build with OpenMesh support" )
      include_directories( SYSTEM ${OPENMESH_INCLUDE_DIRS} )
      list( APPEND SERVICE_LIBS ${OPENMESH_LIBRARIES} )
   else()
      set( WALBERLA_BUILD_WITH_OPENMESH OFF CACHE BOOL "Build with OpenMesh support" FORCE )
   endif()
endif()

#make python available
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

#make eigen SuiteSparse available
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR};${CAMKE_MODULE_PATH}")
find_package(SuiteSparse)
include_directories(${SuiteSparse_INCLUDE_DIRS})

aux_source_directory(. src)
file(GLOB * header)

link_directories("/usr/lib/x86_64-linux-gnu/")

add_executable(data_aug ${src} ${header})
target_link_libraries(data_aug GLU glut GL ${Boost_LIBRARIES} ${SuiteSparse_LIBRARIES} ${ASSIMP_LIBRARIES} ${OPENMESH_LIBRARIES} ${PYTHON_LIBRARIES})

