cmake_minimum_required(VERSION 2.6)
project(SpLLT)
enable_language(Fortran)
enable_language(C)

#########################################
# Version number
set(SPLLT_VERSION_MAJOR 1)
set(SPLLT_VERSION_MINOR 0)
set(SPLLT_VERSION_MICRO 0)


# Include our own cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" )

##################################################################################
# Options

#########################################
# Runtime system

set(RUNTIME "STF" CACHE STRING "Runtime system")
set_property(CACHE RUNTIME PROPERTY STRINGS STF StarPU OMP Parsec)
if(${RUNTIME} MATCHES "STF")
  message(STATUS "Build sequential version of SpLLT")
else()
  message(STATUS "Build SpLLT with ${RUNTIME} runtime system")
endif()

#########################################
# Enable GPU
option(USE_GPU "Enable GPU" OFF)
message(STATUS "Enable GPU: " ${USE_GPU})   

#########################################
# Enable MPI
option(USE_MPI "Enable MPI" OFF)
message(STATUS "Enable MPI: " ${USE_MPI})

#########################################
# Build for KNL
option(USE_KNL "Run on KNL device" OFF)
message(STATUS "Run on KNL device: " ${USE_KNL})

# Test driver
## SpLLT: Cholesky solver
## SpLDLT: LDLT solver
# set(TEST_DRIVER "SpLLT" CACHE STRING "Test driver")
# set_property(CACHE TEST_DRIVER PROPERTY STRINGS SpLLT SpLDLT)
# MESSAGE( STATUS "Test driver:         " ${TEST_DRIVER} )

# MA87 test driver
option(MA87_TEST_DRIVER "Use the ma87 test driver" OFF)
message(STATUS "Use MA87 test driver: " ${MA87_TEST_DRIVER})

##################################################################################
# Compiler option

#########################################
# Custom compiler flags

# Set C flags
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2" )
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2 -ffpe-trap=underflow,denormal" )

# Set CXX flags
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2" )
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -ffpe-trap=underflow,denormal" )

# Set Fortran flags
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O2" )
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O2 -ffpe-trap=zero,underflow,denormal")

#########################################
# Build type

# Default build type
set(default_build_type "RelWithDebInfo")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Using default build type '${default_build_type}' because none was specified")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build (Debug, Release, MinSizeRel, RelWithDebInfo)" FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_Fortran_FLAGS_ALL ${CMAKE_Fortran_FLAGS})
set(CMAKE_C_FLAGS_ALL ${CMAKE_C_FLAGS})

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  list(APPEND CMAKE_Fortran_FLAGS_ALL ${CMAKE_Fortran_FLAGS_DEBUG})
  list(APPEND CMAKE_C_FLAGS_ALL ${CMAKE_C_FLAGS_DEBUG})
elseif(CMAKE_BUILD_TYPE MATCHES "Release")
  list(APPEND CMAKE_Fortran_FLAGS_ALL ${CMAKE_Fortran_FLAGS_RELEASE})
  list(APPEND CMAKE_C_FLAGS_ALL ${CMAKE_C_FLAGS_RELEASE})
elseif(CMAKE_BUILD_TYPE MATCHES "MinSizeRel")
  list(APPEND CMAKE_Fortran_FLAGS_ALL ${CMAKE_Fortran_FLAGS_MINSIZEREL})
  list(APPEND CMAKE_C_FLAGS_ALL ${CMAKE_C_FLAGS_MINSIZEREL})
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")  
  list(APPEND CMAKE_Fortran_FLAGS_ALL ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO})
  list(APPEND CMAKE_C_FLAGS_ALL ${CMAKE_C_FLAGS_RELWITHDEBINFO})
endif()

# message( STATUS "Fortran FLAGS: " ${CMAKE_Fortran_FLAGS} )


##################################################################################
# Third-party libraries

#########################################
# BLAS
set(BLAS_LAPACK "MKL" CACHE STRING "BLAS and LAPACK libraries")
set_property(CACHE BLAS_LAPACK PROPERTY STRINGS Ref MKL)
message( STATUS "BLAS/LAPACK: " ${BLAS_LAPACK} )

#########################################
# SPRAL
set(SPRAL_DIR $ENV{SPRALDIR} CACHE PATH "Location of Spral Library")
include_directories(${SPRAL_DIR})
message(STATUS "SPRAL Directory: " ${SPRAL_DIR})   

#########################################
# Runtime system

# add spLLT test driver
if (${RUNTIME} MATCHES "STF")

  # Use STF code
  add_definitions(-DSPLLT_USE_STF)

elseif(${RUNTIME} MATCHES "StarPU")

  # Use nested STF
  option(SPLLT_USE_NESTED_STF "Use nested STF" OFF)

  if (SPLLT_USE_NESTED_STF MATCHES ON)
    # nested STF
    add_definitions(-DSPLLT_USE_NESTED_STF)
  endif()

  # message( STATUS "Nested STF:      " ${SPLLT_USE_NESTED_STF} )

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")

  # Use StarPU code
  add_definitions(-DSPLLT_USE_STARPU)

  # StarPU lib
  set(STARPU_DIR $ENV{STARPUDIR} CACHE PATH "Location of StarPU directory")
  set(STARPU_INCLUDE ${STARPU_DIR}/include/starpu/1.3 CACHE PATH "Location of StarPU includes")
  include_directories(${STARPU_INCLUDE})
  find_library(STARPU_LIBRARY starpu-1.3 PATHS ${STARPU_DIR}/lib)
  
  #########################################
  # HWLOC
  set(HWLOC_DIR $ENV{HWLOCDIR} CACHE PATH "Location of HWLOC Library")
  include_directories(${HWLOC_DIR}/include)  
  find_library(HWLOC_LIBRARY hwloc PATHS ${HWLOC_DIR}/lib)

elseif(${RUNTIME} MATCHES "OMP")

  # TODO find_package omp
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")

  # Enable OMP
  add_definitions(-DSPLLT_USE_OMP)

elseif(${RUNTIME} MATCHES "Parsec")

  # Set C flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -mcx16" )
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 -mcx16" )

  # Parsec
  set(PARSEC_DIR $ENV{PARSECDIR} CACHE PATH "Location of Parsec Library")
  set(PARSEC_SRC_DIR $ENV{PARSECSRCDIR} CACHE PATH "Location of Parsec sources")
  include_directories(${PARSEC_SRC_DIR})
  include_directories(${PARSEC_SRC_DIR}/parsec)
  include_directories(${PARSEC_SRC_DIR}/parsec/include)
  include_directories(${PARSEC_SRC_DIR}/parsec/include/dague)
  include_directories(${PARSEC_DIR}/parsec/include)
  include_directories(${PARSEC_DIR}/parsec/fortran)
  include_directories(${PARSEC_DIR}/data_dist/matrix)

  message(STATUS "Parsec directory ENV: " $ENV{PARSECDIR})
  message(STATUS "Parsec directory: " ${PARSEC_SRC_DIR})
  
  # find_library(PARSEC_LIBRARY dague_profilef PATHS ${PARSEC_DIR}/dague/fortran)
  # find_library(PARSEC_LIBRARY daguef PATHS ${PARSEC_DIR}/dague/fortran)
  # find_library(PARSEC_LIBRARY dague_distribution PATHS ${PARSEC_DIR}/data_dist)
  # find_library(PARSEC_LIBRARY dague PATHS ${PARSEC_DIR}/dague)
  # find_library(PARSEC_LIBRARY dague-base PATHS ${PARSEC_DIR}/dague)

  #########################################
  # Parsec libraries
  find_library(PARSEC_FORTRAN_LIBRARY
    NAMES parsecf
    PATHS ${PARSEC_DIR}/parsec/fortran)

  find_library(PARSEC_PROFILEF_LIBRARY
    NAMES parsec_profilef 
    PATHS ${PARSEC_DIR}/parsec/fortran)

  find_library(PARSEC_DIST_LIBRARY
    NAMES parsec_distribution
    PATHS ${PARSEC_DIR}/data_dist)

  find_library(PARSEC_DIST_MAT_LIBRARY
    NAMES parsec_distribution_matrix
    PATHS ${PARSEC_DIR}/data_dist/matrix)

  find_library(PARSEC_BASE_LIBRARY
    NAMES parsec-base
    PATHS ${PARSEC_DIR}/parsec)

  find_library(PARSEC_LIB
    NAMES parsec
    PATHS ${PARSEC_DIR}/parsec)

  set (PARSEC_LIBRARY ${PARSEC_PROFILEF_LIBRARY}
    ${PARSEC_FORTRAN_LIBRARY} ${PARSEC_DIST_MAT_LIBRARY} 
    ${PARSEC_DIST_LIBRARY} ${PARSEC_LIB} ${PARSEC_BASE_LIBRARY})

  # message(STATUS "Parsec library: ${PARSEC_LIBRARY}")

  # Dague PP
  # set(DAGUEPP_CFLAGS "--noline" CACHE STRING "Additional daguepp precompiling flags" )
  # set(PARSECPP $ENV{DAGUEPP} CACHE PATH "Parsec JDF compiler") 
  # set(daguepp_EXE ${PARSECPP})

  # Find and configure JDF compiler parsec_ptgpp 
  set(PARSEC_PTGPP_CFLAGS "--noline" CACHE STRING "Additional daguepp precompiling flags")
  set(PARSEC_PTGPP $ENV{PARSECPP} CACHE PATH "Parsec JDF compiler") 
  set(parsec_ptgpp_EXE ${PARSEC_PTGPP})
  
  # message(STATUS "Dague PP: ${daguepp_EXE}")

  #########################################
  # HWLOC
  set(HWLOC_DIR $ENV{HWLOCDIR} CACHE PATH "Location of HWLOC Library")
  include_directories(${HWLOC_DIR}/include)  
  find_library(HWLOC_LIBRARY hwloc PATHS ${HWLOC_DIR}/lib)

  #########################################
  # PAPI
  find_library(PAPI_LIBRARY papi)
  message(STATUS "PAPI: " ${PAPI_LIBRARY})

  # Enable Parsec
  add_definitions(-DSPLLT_USE_PARSEC)

endif()

if(${USE_GPU} MATCHES ON)

  # Use GPU
  add_definitions(-DSPLLT_USE_GPU)

  # set(CUDA_NVCC_FLAGS "-g -arch=compute_20 -code=compute_20,sm_20,sm_35" CACHE
  #   STRING "NVCC flags")

  #########################################
  # CUDA
  find_package(CUDA REQUIRED)

  # set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-g;-arch=compute_20,code=compute_20,sm_20,sm_35")
  # set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-g;-arch=compute_20;-code=compute_20,sm_20,sm_35" CACHE STRING "NVCC flags")

  # set(CUDA_NVCC_FLAGS "TET" CACHE STRING "NVCC flags")
  list(APPEND CUDA_NVCC_FLAGS "-g")
  
  MESSAGE( STATUS "CUDA found:         " ${CUDA_FOUND} )
  MESSAGE( STATUS "CUDA host compiler: " ${CUDA_HOST_COMPILER} )
  MESSAGE( STATUS "CUDA nvcc compiler: " ${CUDA_NVCC_EXECUTABLE})
  MESSAGE( STATUS "NVCC FLAGS:         " ${CUDA_NVCC_FLAGS} )

  # enable_language(CUDA)
  # include(FindCUDA)

  # Allow for dynamic parallelism with CUBLAS
  # set(CUDA_SEPARABLE_COMPILATION TRUE)

  # MAGMA library
  set(MAGMA_DIR $ENV{MAGMADIR} CACHE PATH "Location of MAGMA Library")
  set(MAGMA_SRCDIR $ENV{MAGMASRCDIR} CACHE PATH "Location of MAGMA Sources")
  include_directories(${MAGMA_DIR}/include)
  include_directories(${MAGMA_SRCDIR}/control)     
  find_library(MAGMA_LIBRARY magma PATHS ${MAGMA_DIR}/lib)
endif()

#########################################
# MPI
if(${USE_MPI} MATCHES ON)

  # Use GPU
  add_definitions(-DSPLLT_USE_MPI)

  set(MPI_DIR $ENV{MPIDIR} CACHE PATH "Location of MPI Library")

  include_directories(${MPI_DIR}/include)
  find_library(MPI_LIBRARY mpi PATHS ${MPI_DIR}/lib)

endif()


# Set flags for the MA87 test driver

if(MA87_TEST_DRIVER MATCHES ON)
  # TODO find_package omp
  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
endif()

# Adds sources files

add_subdirectory(src)

# Generate the MA87 test driver

if(MA87_TEST_DRIVER MATCHES ON)

  MESSAGE( STATUS "MA87 test driver:         " ${MA87_TEST_DRIVER} )

  add_executable(spllt_test drivers/run_ma87.f90)
else()
  # add_executable(spllt_test spllt_test.F90)
  add_executable(spllt_test drivers/spllt_test.F90)
endif()

target_include_directories(spllt_test PUBLIC ${CMAKE_BINARY_DIR}/src)
target_link_libraries(spllt_test spllt)
# find libraries

if(${RUNTIME} MATCHES "StarPU")
  # StarPU
  target_link_libraries(spllt_test ${STARPU_LIBRARY})
  # HWLOC
  target_link_libraries(spllt_test ${HWLOC_LIBRARY})

elseif(${RUNTIME} MATCHES "Parsec")

  # Parsec
  target_link_libraries(spllt_test ${PARSEC_LIBRARY})

  # HWLOC
  target_link_libraries(spllt_test ${HWLOC_LIBRARY})

  # PAPI
  target_link_libraries(spllt_test ${PAPI_LIBRARY})

endif()

## GPU
if(${USE_GPU} MATCHES ON)
  # MAGMA
  target_link_libraries(spllt_test ${MAGMA_LIBRARY})

  target_link_libraries(spllt_test spllt_cuda)
endif()

## MPI
if(${USE_MPI} MATCHES ON)
  target_link_libraries(spllt_test ${MPI_LIBRARY})
endif()

# spral
target_include_directories(spllt_test PUBLIC ${SPRAL_DIR})
find_library(SPRAL_LIBRARY spral PATHS ${SPRAL_DIR})
target_link_libraries(spllt_test ${SPRAL_LIBRARY})
# MESSAGE( STATUS "spral dir:         " ${SPRAL_DIR})
# MESSAGE( STATUS "spral lib:         " ${SPRAL_LIBRARY})

#########################################
# Metis
set(METIS_DIR $ENV{METISDIR} CACHE PATH "Location of Metis Library")
find_library(METIS_LIBRARY metis PATHS ${METIS_DIR})
target_link_libraries(spllt_test ${METIS_LIBRARY})

# HWLOC
set(HWLOC_DIR $ENV{HWLOCDIR} CACHE PATH "Location of HWLOC Library")
include_directories(${HWLOC_DIR}/include)  
find_library(HWLOC_LIBRARY hwloc PATHS ${HWLOC_DIR}/lib)
target_link_libraries(spllt_test ${HWLOC_LIBRARY})

# BLAS and LAPACK
if (${BLAS_LAPACK} MATCHES "MKL")

  # MKL library

  if(${USE_KNL} MATCHES ON)

    find_library(MKL_IFACE_LIBRARY mkl_intel_lp64 PATHS ${MKL_LIBS} ${MKL_DIR}/lib/intel64)
    find_library(MKL_THREAD_LIBRARY mkl_intel_thread PATHS ${MKL_ROOT}/lib/intel64)
    
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")

    target_link_libraries(spllt_test iomp5)

  else()

    set(MKL_LIBS $ENV{MKL_LIBS})
    set(MKL_ROOT "/opt/intel/mkl" CACHE PATH "Location of MKL Library")
    set(MKL_DIR $ENV{MKLROOT})

    find_library(MKL_IFACE_LIBRARY mkl_gf_lp64 PATHS ${MKL_LIBS} ${MKL_ROOT}/lib/intel64 ${MKL_DIR}/lib/intel64)
    # find_library(MKL_THREAD_LIBRARY mkl_gnu_thread PATHS ${MKL_ROOT}/lib/intel64)
    find_library(MKL_SEQ_LIBRARY mkl_sequential PATHS ${MKL_LIBS} ${MKL_DEFAULT_ROOT}/lib/intel64 ${MKL_DIR}/lib/intel64)
    find_library(MKL_CORE_LIBRARY mkl_core PATHS ${MKL_LIBS} ${MKL_ROOT}/lib/intel64 ${MKL_DIR}/lib/intel64) 
    # set(MKL_LIBRARIES ${MKL_IFACE_LIBRARY} ${MKL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY})

    MESSAGE( STATUS "MKL Root: " ${MKL_DIR} )
    MESSAGE( STATUS "MKL MKL_IFACE_LIBRARY: " ${MKL_IFACE_LIBRARY} )

  endif()

  set(BLAS_LAPACK_LIBRARIES ${MKL_IFACE_LIBRARY} ${MKL_SEQ_LIBRARY} ${MKL_CORE_LIBRARY})
  MESSAGE( STATUS "BLAS and LAPACK Libraries: " ${BLAS_LAPACK_LIBRARIES} )

else()

  # Reference BLAS libraries

  find_library(BLAS_LIBRARY blas)
  find_library(LAPACK_LIBRARY lapack)
  set(BLAS_LAPACK_LIBRARIES ${BLAS_LIBRARY} ${LAPACK_LIBRARY})

endif()

target_link_libraries(spllt_test ${BLAS_LAPACK_LIBRARIES})

# pthread
target_link_libraries(spllt_test pthread)

## std C++ (needed for linking with ssids)
target_link_libraries(spllt_test stdc++)

target_link_libraries(spllt_test rt)

#########################################
# Print summary

message("Configuration done\n"
  "Build type: ${CMAKE_BUILD_TYPE}\n"
  "Fortran compiler: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})\n"
  "Fortran compiler flags: ${CMAKE_Fortran_FLAGS_ALL}\n"
  "C compiler:  ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})\n"
  "C compiler flags: ${CMAKE_Fortran_FLAGS_ALL}\n")
