cmake_minimum_required(VERSION 2.4)
project(UnigenQA)
set(CMAKE_cpp_STANDARD 11)
set(CMAKE_BUILD_TYPE debug)

set(EXE_SOURCES
        RunUnigenQA.C)

set(UNIGEN_SOURCES
        URun.cxx
        UEvent.cxx
        UParticle.cxx)

set(QA_SOURCES
        UnigenQA.cxx)

string(REPLACE ".cxx" ".h" UNIGEN_HEADERS "${UNIGEN_SOURCES}")
string(REPLACE ".cxx" ".h" QA_HEADERS "${QA_SOURCES}")

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Telling CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/lib")

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Core MathCore MathMore RIO Hist Tree Net TreePlayer)

# Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
add_definitions(${ROOT_CXX_FLAGS})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})

ROOT_GENERATE_DICTIONARY(G__Unigen ${UNIGEN_HEADERS} LINKDEF UnigenLinkDef.h)
add_library(Unigen SHARED ${UNIGEN_SOURCES} G__Unigen)

add_library(UnigenQA SHARED ${QA_SOURCES})
target_link_libraries(UnigenQA ${ROOT_LIBRARIES} Unigen)

# Create a main program using the library
add_executable(RunUnigenQA ${EXE_SOURCES})
target_link_libraries(RunUnigenQA UnigenQA ${ROOT_LIBRARIES})
target_include_directories(RunUnigenQA PUBLIC)

