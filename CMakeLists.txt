
cmake_minimum_required(VERSION 3.8)

option(CXX "enable C++ compilation" ON)
if (CXX)
  enable_language(CXX)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILE_FEATURES cxx_std_17)

project(PRactIP CXX C)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
add_definitions(-DPACKAGE_VERSION="0.0.2")

# find_package(PkgConfig REQUIRED)

# pkg_check_modules(VIENNARNA REQUIRED RNAlib2)
# include_directories(${VIENNARNA_INCLUDE_DIRS})
# add_definitions(-DHAVE_VIENNA18)
# add_definitions(-DHAVE_VIENNA20)
find_library(PTHREAD_LIBRARY pthread)

if(ENABLE_CPLEX)
  find_package(CPLEX)
  if (CPLEX_FOUND)
    add_definitions(-DWITH_CPLEX)
    include_directories(${CPLEX_INCLUDE_DIR})
  endif()
elseif(ENABLE_GUROBI)
  find_package(GUROBI)
  if (GUROBI_FOUND)
    add_definitions(-DWITH_GUROBI)
    include_directories(${GUROBI_INCLUDE_DIR})
  endif()
else()
  find_package(GLPK 4.33)
  if(GLPK_FOUND)
    find_library(ZLIB_LIBRARY z)
    find_library(GMP_LIBRARY gmp)
    find_library(DL_LIBRARY dl)
    find_library(LTDL_LIBRARY ltdl)
    find_library(MPFR_LIBRARY mpfr)
    add_definitions(-DWITH_GLPK)
    include_directories(${GLPK_INCLUDE_DIR})
  endif()
endif()

add_executable(
  practip
  src/practip.cpp
  src/practip.h
  src/ip.cpp
  src/ip.h
  src/typedef.h
  #src/defparam.dat
  src/cxxopts.hpp
  )
target_link_libraries(practip 
  ${VIENNARNA_LDFLAGS} ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${MPFR_LIBRARY}
  ${CPLEX_LIBRARY} ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} ${PTHREAD_LIBRARY} ${CMAKE_DL_LIBS})
target_compile_features(practip PUBLIC cxx_std_17)

install(TARGETS practip RUNTIME DESTINATION bin)
