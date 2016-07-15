cmake_minimum_required (VERSION 2.6)
project (SpLib C CXX Fortran)

# Include our own cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set C flags 
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -g -O2")

# Set Fortran flags 
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O2")

# add_library(hsl_ma87 hsl_ma87d.F90 
#   $ENV{HSLPACKDIR}/hsl_mc78/hsl_mc78i.f90
#   $ENV{HSLPACKDIR}/hsl_mc34/hsl_mc34d.f90)

SET(RUNTIME "STF" CACHE STRING "Runtime system")
set_property(CACHE RUNTIME PROPERTY STRINGS STF StarPU OMP Parsec)
MESSAGE( STATUS "Runtime:         " ${RUNTIME} )

# Use GPU
option(USE_GPU "Use GPU" OFF)
MESSAGE( STATUS "Use GPU:         " ${USE_GPU} )   

# Use MPI
option(USE_MPI "Use MPI" OFF)
MESSAGE( STATUS "Use MPI:         " ${USE_MPI} )   

# add spLLT test driver

if (${RUNTIME} MATCHES "STF")

  # Use STF code
  add_definitions(-DSPLLT_USE_STF)

elseif(${RUNTIME} MATCHES "StarPU")

  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")

  # Use StarPU code
  add_definitions(-DSPLLT_USE_STARPU)

  # StarPU lib
  set(STARPU_DIR $ENV{STARPUDIR} CACHE PATH "Location of StarPU Library")
  include_directories(${STARPU_DIR}/include/starpu/1.3)
  find_library(STARPU_LIBRARY starpu-1.3 PATHS ${STARPU_DIR}/lib)
  
  # HWLOC
  set(HWLOC_DIR $ENV{HWLOCDIR} CACHE PATH "Location of HWLOC Library")
  include_directories(${HWLOC_DIR}/include)  
  find_library(HWLOC_LIBRARY hwloc PATHS ${HWLOC_DIR}/lib)

elseif(${RUNTIME} MATCHES "OMP")

  set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")

  # Use OMP code
  add_definitions(-DSPLLT_USE_OMP)

endif()

if(${USE_GPU} MATCHES ON)

  # Use GPU
  add_definitions(-DSPLLT_USE_GPU)

  # set(CUDA_NVCC_FLAGS "-g -arch=compute_20 -code=compute_20,sm_20,sm_35" CACHE
  #   STRING "NVCC flags")

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

add_subdirectory(src)

add_executable(spllt_test spllt_test.F90)
target_include_directories(spllt_test PUBLIC ${CMAKE_BINARY_DIR}/src)
target_link_libraries(spllt_test spllt)

# find libraries

if(${RUNTIME} MATCHES "StarPU")
  # StarPU
  target_link_libraries(spllt_test ${STARPU_LIBRARY})
  # HWLOC
  target_link_libraries(spllt_test ${HWLOC_LIBRARY})
endif()

if(${USE_GPU} MATCHES ON)
  # MAGMA
  target_link_libraries(spllt_test ${MAGMA_LIBRARY})

  target_link_libraries(spllt_test spllt_cuda)
endif()

# Metis
set(METIS_DIR $ENV{METISDIR} CACHE PATH "Location of Metis Library")
find_library(METIS_LIBRARY metis PATHS ${METIS_DIR})
target_link_libraries(spllt_test ${METIS_LIBRARY})

# spral
set(SPRAL_DIR $ENV{SPRALDIR} CACHE PATH "Location of Spral Library")
target_include_directories(spllt_test PUBLIC ${SPRAL_DIR})
find_library(SPRAL_LIBRARY spral PATHS ${SPRAL_DIR} ${SPRAL_DIR}/lib)
target_link_libraries(spllt_test ${SPRAL_LIBRARY})

# BLAS and LAPACK
set(MKL_LIBS $ENV{MKL_LIBS})
set(MKL_ROOT "/opt/intel/mkl" CACHE PATH "Location of MKL Library")
set(MKL_DIR $ENV{MKLROOT})

find_library(MKL_IFACE_LIBRARY mkl_gf_lp64 PATHS ${MKL_LIBS} ${MKL_ROOT}/lib/intel64 ${MKL_DIR}/lib/intel64)
# find_library(MKL_THREAD_LIBRARY mkl_gnu_thread PATHS ${MKL_ROOT}/lib/intel64)
find_library(MKL_SEQ_LIBRARY mkl_sequential PATHS ${MKL_LIBS} ${MKL_DEFAULT_ROOT}/lib/intel64 ${MKL_DIR}/lib/intel64)
find_library(MKL_CORE_LIBRARY mkl_core PATHS ${MKL_LIBS} ${MKL_ROOT}/lib/intel64 ${MKL_DIR}/lib/intel64) 
# set(MKL_LIBRARIES ${MKL_IFACE_LIBRARY} ${MKL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY})
set(MKL_LIBRARIES ${MKL_IFACE_LIBRARY} ${MKL_SEQ_LIBRARY} ${MKL_CORE_LIBRARY})

target_link_libraries(spllt_test ${MKL_LIBRARIES})

# pthread
target_link_libraries(spllt_test pthread)
