cmake_minimum_required(VERSION 2.6)
include(CheckFunctionExists)
include(CheckFortranFunctionExists)
project(SpLLT)
enable_language(Fortran)
enable_language(C)

########################################
# Version number
set(SPLLT_VERSION_MAJOR 1)
set(SPLLT_VERSION_MINOR 0)
set(SPLLT_VERSION_MICRO 0)

# Include our own cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" )

################################################################################
# Options

########################################
# Runtime system

set(RUNTIME "STF" CACHE STRING "Runtime system")
set_property(CACHE RUNTIME PROPERTY STRINGS STF StarPU OMP Parsec)
if(${RUNTIME} MATCHES "STF")
  message(STATUS "Build SEQUENTIAL version of SpLLT")
else()
  message(STATUS "Build SpLLT with ${RUNTIME} runtime system")
endif()

########################################
# Enable GPU
option(USE_GPU "Enable GPU" OFF)
message(STATUS "Enable GPU: " ${USE_GPU})   

########################################
# Enable MPI
option(USE_MPI "Enable MPI" OFF)
message(STATUS "Enable MPI: " ${USE_MPI})

########################################
# Build for KNL
option(USE_KNL "Run on KNL device" OFF)
message(STATUS "Run on KNL device: " ${USE_KNL})

# Test driver
## SpLLT: Cholesky solver
## SpLDLT: LDLT solver
# set(TEST_DRIVER "SpLLT" CACHE STRING "Test driver")
# set_property(CACHE TEST_DRIVER PROPERTY STRINGS SpLLT SpLDLT)
# MESSAGE( STATUS "Test driver:         " ${TEST_DRIVER} )

################################################################################
# Compiler options

########################################
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

########################################
# Build type

# Set default build type
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
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all -Wall")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -fcheck=all -Wall")
# CACHE STRING "Flags used by the compiler during release builds with debug info." FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -mtune=native")
# CACHE STRING "Flags used by the compiler during release builds." FORCE)

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

################################################################################
# Third-party libraries

set(LIBS)

########################################
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
  
  ########################################
  # StarPU
  find_package(STARPU)
  if (STARPU_FOUND)
    include_directories( ${STARPU_INCLUDE_DIRS} )
    link_directories( ${STARPU_LIBRARY_DIRS} )
    link_libraries( ${STARPU_LIBRARIES} )
    set(LIBS ${LIBS} ${STARPU_LIBRARIES})
    find_path(STARPU_F_MOD_FOUND fstarpu_mod.f90
              HINTS ${STARPU_INCLUDE_DIRS}
             )
    set(STARPU_F_MOD "${STARPU_F_MOD_FOUND}/fstarpu_mod.f90")
    include_directories(${STARPU_INCLUDE_DIRS}) 
  else()
    message(FATAL_ERROR "StarPU not found!")
  endif(STARPU_FOUND)
  
  # Use StarPU code
  add_definitions(-DSPLLT_USE_STARPU)

  ########################################
  # HWLOC
  set(HWLOC_DIR $ENV{HWLOCDIR} CACHE PATH "Location of HWLOC Library")
  include_directories(${HWLOC_DIR}/include)  
  find_library(HWLOC_LIBRARY hwloc PATHS ${HWLOC_DIR}/lib)
  set(LIBS ${LIBS} ${HWLOC_LIBRARIES})

elseif(${RUNTIME} MATCHES "OMP")
  include_directories("./src/include/")

  # TODO find_package omp
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
  else()
    set(ERROR " Unrecognised compiler type so that no openmp flag used")
  endif()

  # Enable OMP
  add_definitions(-DSPLLT_USE_OMP)

  option(WITH_SKERNEL_SCATTER "Fast timer on scatter step in solve kernels" OFF)
  if(${WITH_SKERNEL_SCATTER} MATCHES ON)
    add_definitions(-DSPLLT_SOLVE_KERNEL_SCATTER)
    message(STATUS "WITH_SKERNEL_SCATTER ... yes")
  else()
    message(STATUS "WITH_SKERNEL_SCATTER ... no")
  endif()
  
  option(WITH_SKERNEL_GATHER "Fast timer on gather step in solve kernels" OFF)
  if(${WITH_SKERNEL_GATHER} MATCHES ON)
    add_definitions(-DSPLLT_SOLVE_KERNEL_GATHER)
    message(STATUS "WITH_SKERNEL_GATHER ... yes")
  else()
    message(STATUS "WITH_SKERNEL_GATHER ... no")
  endif()
  
  option(WITH_PROFILING "Activate the flop count" OFF)
  if(${WITH_PROFILING} MATCHES ON)
    add_definitions(-DSPLLT_PROFILING_FLOP)
    message(STATUS "WITH_PROFILING ... yes")
  else()
    message(STATUS "WITH_PROFILING ... no")
  endif()
  
  option(WITH_TRACE "Activate the trace" OFF)
  if(${WITH_TRACE} MATCHES ON)
    add_definitions(-DSPLLT_OMP_TRACE)
    message(STATUS "WITH_TRACE for OMP ... yes")
  else()
    message(STATUS "WITH_TRACE for OMP ... no")
  endif()
  
  option(DRIVER_CHKERR "Check the backward error" OFF)
  if(${DRIVER_CHKERR} MATCHES ON)
    add_definitions(-DSPLLT_DRIVER_CHECK_ERROR)
    message(STATUS "DRIVER_CHKERR ... yes")
  else()
    message(STATUS "DRIVER_CHKERR ... no")
  endif()
  
  option(VERBOSE_TASK "Display the dependencies of the task during the execution" OFF)
  if(${VERBOSE_TASK} MATCHES ON)
    add_definitions(-DSPLLT_VERBOSE)
    message(STATUS "VERBOSE_TASK ... yes")
  else()
    message(STATUS "VERBOSE_TASK ... no")
  endif()
  
  option(WITH_TIMER_TASK "Record the time spend in tasks" OFF)
  if(${WITH_TIMER_TASK} MATCHES ON)
    add_definitions(-DSPLLT_TIMER_TASKS)
    message(STATUS "WITH_TIMER_TASK ... yes")
  else()
    message(STATUS "WITH_TIMER_TASK ... no")
  endif()
  
  option(WITH_TIMER_TASK_SUBMIT "Record the time spend to submit tasks" OFF)
  if(${WITH_TIMER_TASK_SUBMIT} MATCHES ON)
    add_definitions(-DSPLLT_TIMER_TASKS_SUBMISSION)
    message(STATUS "WITH_TIMER_TASK_SUBMIT ... yes")
  else()
    message(STATUS "WITH_TIMER_TASK_SUBMIT ... no")
  endif()

  option(WITH_TASK_LOCK "Activate a lock that ensure the execution of the
  solve tasks is nperformed before the submission of the last task" OFF)
  if(${WITH_TASK_LOCK} MATCHES ON)
    add_definitions(-DSOLVE_TASK_LOCKED)
    message(STATUS "WITH_TASK_LOCK ... yes")
  else()
    message(STATUS "WITH_TASK_LOCK ... no")
  endif()

elseif(${RUNTIME} MATCHES "Parsec")

  # Set C flags
  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -mcx16" )
  # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 -mcx16" )

  ########################################
  # Parsec
  set(PARSEC_DIR $ENV{PARSEC_DIR} CACHE PATH "Location of Parsec Library")
  set(PARSEC_SRC_DIR $ENV{PARSEC_SRCDIR} CACHE PATH "Location of Parsec sources")
  include_directories(${PARSEC_SRC_DIR})
  include_directories(${PARSEC_SRC_DIR}/parsec)
  include_directories(${PARSEC_SRC_DIR}/parsec/include)
  include_directories(${PARSEC_SRC_DIR}/parsec/include/parsec)
  include_directories(${PARSEC_DIR}/parsec/include)
  include_directories(${PARSEC_DIR}/parsec/mod_files)
  include_directories(${PARSEC_DIR}/data_dist/matrix)

  # message(STATUS "Parsec directory ENV: " $ENV{PARSECDIR})
  message(STATUS "Parsec source directory: " ${PARSEC_SRC_DIR})
  
  ########################################
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

  set (PARSEC_LIBRARIES ${PARSEC_PROFILEF_LIBRARY}
    ${PARSEC_FORTRAN_LIBRARY} ${PARSEC_DIST_MAT_LIBRARY} 
    ${PARSEC_DIST_LIBRARY} ${PARSEC_LIB} ${PARSEC_BASE_LIBRARY})
  set(LIBS ${LIBS} ${PARSEC_LIBRARIES})

  # message(STATUS "Parsec library: ${PARSEC_LIBRARY}")

  # Dague PP
  # set(DAGUEPP_CFLAGS "--noline" CACHE STRING "Additional daguepp precompiling flags" )
  # set(PARSECPP $ENV{DAGUEPP} CACHE PATH "Parsec JDF compiler") 
  # set(daguepp_EXE ${PARSECPP})

  # Find and configure JDF compiler parsec_ptgpp 
  set(PARSEC_PTGPP_CFLAGS "--noline" CACHE STRING "Additional daguepp precompiling flags")
  set(PARSEC_PTGPP $ENV{PARSECPP} CACHE PATH "Parsec JDF compiler") 
  set(parsec_ptgpp_EXE ${PARSEC_PTGPP})
  
  # # message(STATUS "Dague PP: ${daguepp_EXE}")

  # find_package(PARSEC)
  # if (PARSEC_FOUND)
  #   include_directories(${PARSEC_INCLUDE_DIRS})
  #   link_directories(${PARSEC_LIBRARY_DIRS})
  #   link_libraries(${PARSEC_LIBRARIES})
  #   # Adds directory containing module files
  #   include_directories(${PARSEC_LIBRARY_DIRS}/mod_files)
  #   message("Parsec inc dirs ${PARSEC_INCLUDE_DIRS}")
  #   message("Parsec lib dirs ${PARSEC_LIBRARY_DIRS}")
  #   set(LIBS ${LIBS} ${PARSEC_LIBRARIES})
  #   # Adds JDF compiler and flags
  #   set(PARSEC_PTGPP_CFLAGS "--noline" CACHE STRING "Additional ptgpp precompiling flags")
  #   set(parsec_ptgpp_EXE ${PARSEC_PARSEC_PTGPP})
  # else()
  #   message(FATAL_ERROR "PARSEC not found!")
  # endif(PARSEC_FOUND)

  ########################################
  # HWLOC
  set(HWLOC_DIR $ENV{HWLOCDIR} CACHE PATH "Location of HWLOC Library")
  include_directories(${HWLOC_DIR}/include)  
  find_library(HWLOC_LIBRARY hwloc PATHS ${HWLOC_DIR}/lib)
  set(LIBS ${LIBS} ${HWLOC_LIBRARIES})

  ########################################
  # PAPI
  find_library(PAPI_LIBRARY papi)
  # message(STATUS "PAPI: " ${PAPI_LIBRARY})
  set(LIBS ${LIBS} ${PAPI_LIBRARY})

  # Enable Parsec
  add_definitions(-DSPLLT_USE_PARSEC)

endif()

if(${USE_GPU} MATCHES ON)

  # Use GPU
  add_definitions(-DSPLLT_USE_GPU)

  # set(CUDA_NVCC_FLAGS "-g -arch=compute_20 -code=compute_20,sm_20,sm_35" CACHE
  #   STRING "NVCC flags")

  ########################################
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
  find_library(MAGMA_LIBRARIES magma PATHS ${MAGMA_DIR}/lib)

  set(LIBS ${LIBS} ${MAGMA_LIBRARIES})

  set(LIBS ${LIBS} spllt_cuda)

endif()

########################################
# MPI
if(${USE_MPI} MATCHES ON)

  # Use GPU
  add_definitions(-DSPLLT_USE_MPI)

  set(MPI_DIR $ENV{MPIDIR} CACHE PATH "Location of MPI Library")

  include_directories(${MPI_DIR}/include)
  find_library(MPI_LIBRARY mpi PATHS ${MPI_DIR}/lib)

  set(LIBS ${LIBS} ${PAPI_LIBRARY})

endif()


########################################
# SpLLT
if(${WITH_TRACE} MATCHES ON)
  set(LIBS ${LIBS} spllt_trace)
else()
  set(LIBS ${LIBS} spllt)
endif()

########################################
# SPRAL

# Manage env
if(DEFINED ENV{SPRAL_LIB})
  set(SPRAL_LIB $ENV{SPRAL_LIB} CACHE PATH "Location of SPRAL library")
endif()
if(DEFINED ENV{SPRAL_INC})
  set(SPRAL_INC $ENV{SPRAL_INC} CACHE PATH "Location of SPRAL include")
endif()
if(DEFINED ENV{SPRAL_DIR})
  set(SPRAL_DIR $ENV{SPRAL_DIR} CACHE PATH "Location of SPRAL directory")
endif()

# Include
if(SPRAL_LIB AND SPRAL_INC)
  message(STATUS "SPRAL_INC used : " ${SPRAL_INC})
  message(STATUS "SPRAL_LIB used : " ${SPRAL_LIB})
  include_directories(${SPRAL_INC})
  set(SPRAL_LIB_PATH CACHE PATH ${SPRAL_LIB} "Location of SPRAL library")
else(SPRAL_LIB AND SPRAL_INC)
  if(SPRAL_DIR)
    include_directories(${SPRAL_DIR})
    message(STATUS "SPRAL_DIR used : " ${SPRAL_DIR})
    set(SPRAL_LIB_PATH ${SPRAL_DIR} CACHE PATH "Location of SPRAL library")
  else(SPRAL_DIR)
    message(ERROR " You should provide both SPRAL_LIB and SPRAL_INC, "
      "or SPRAL_DIR.")
  endif(SPRAL_DIR)
endif(SPRAL_LIB AND SPRAL_INC)

# Search library
if(SPRAL_LIB_PATH)
  set(LIB_TO_FIND "spral")

  find_library(SPRAL_LIBRARIES NAMES ${LIB_TO_FIND} PATHS ${SPRAL_LIB_PATH})

  if(SPRAL_LIBRARIES)
    set(LIBS ${LIBS} ${SPRAL_LIBRARIES})
  else(SPRAL_LIBRARIES)
    message(FATAL_ERROR "${SPRAL_LIB_PATH} does not contain lib${LIB_TO_FIND}")
  endif(SPRAL_LIBRARIES)

  unset(LIB_TO_FIND)
else(SPRAL_LIB_PATH)
  message(FATAL_ERROR "No Path for SPRAL library provided. "
    "Define either both SPRAL_LIB and SPRAL_INC, or SPRAL_DIR.")
endif(SPRAL_LIB_PATH)

########################################
# Metis

# Manage env
if(DEFINED ENV{METIS_LIB})
  set(METIS_LIB $ENV{METIS_LIB} CACHE PATH "Location of METIS library")
endif()
if(DEFINED ENV{METIS_INC})
  set(METIS_INC $ENV{METIS_INC} CACHE PATH "Location of METIS include")
endif()
if(DEFINED ENV{METIS_DIR})
  set(METIS_DIR $ENV{METIS_DIR} CACHE PATH "Location of METIS directory")
endif()

# Include
if(METIS_LIB AND METIS_INC)
  message(STATUS "METIS_INC used : " ${METIS_INC})
  message(STATUS "METIS_LIB used : " ${METIS_LIB})
  include_directories(${METIS_INC})
  set(METIS_LIB_PATH CACHE PATH ${METIS_LIB} "Location of METIS library")
else(METIS_LIB AND METIS_INC)
  if(METIS_DIR)
    include_directories(${METIS_DIR}/Lib)
    message(STATUS "METIS_DIR used: " ${METIS_DIR})
    set(METIS_LIB_PATH ${METIS_DIR} CACHE PATH "Location of METIS library")
  # else(METIS_DIR)
  #   message(ERROR " You should provide both METIS_LIB and METIS_INC, "
  #     "or METIS_DIR.")
  endif(METIS_DIR)
endif(METIS_LIB AND METIS_INC)

# Search library

unset(METIS_LIBRARIES)

if(METIS_LIB_PATH)

  find_library(METIS_LIBRARIES NAMES "metis" PATHS ${METIS_LIB_PATH})

else()

  find_library(METIS_LIBRARIES NAMES "metis")

endif(METIS_LIB_PATH)

if(METIS_LIBRARIES)
  set(LIBS ${LIBS} ${METIS_LIBRARIES})
else(METIS_LIBRARIES)
  message(FATAL_ERROR "METIS NOT FOUND")
endif(METIS_LIBRARIES)


#if(METIS_LIB)
#  set(METIS_LIB_PATH CACHE PATH ${METIS_LIB} "Location of METIS library")
#else(METIS_LIB)
#  if(DEFINED ENV{METIS_LIB})
#    set(METIS_LIB_PATH CACHE PATH $ENV{METIS_LIB} "Location of METIS library")
#  else(DEFINED ENV{METIS_LIB})
#    set(METIS_DIR $ENV{METIS_DIR} CACHE PATH "Location of METIS directory")
#    include_directories(${METIS_DIR}/Lib)
#    message(STATUS "METIS Directory: " ${METIS_DIR})
#    set(METIS_LIB_PATH ${METIS_DIR} CACHE PATH "Location of METIS library")
#  endif(DEFINED ENV{METIS_LIB})
#endif(METIS_LIB)
#
#if(METIS_LIB_PATH)
#  set(LIB_TO_FIND "metis")
#
#  find_library(METIS_LIBRARIES NAMES ${LIB_TO_FIND} PATHS ${METIS_LIB_PATH})
#
#  if(METIS_LIBRARIES)
#    set(LIBS ${LIBS} ${METIS_LIBRARIES})
#  else(METIS_LIBRARIES)
#    message(FATAL_ERROR "${METIS_LIB_PATH} does not contain lib${LIB_TO_FIND}")
#  endif(METIS_LIBRARIES)
#
#  unset(LIB_TO_FIND)
#else(METIS_LIB_PATH)
#  message(FATAL_ERROR "No Path for METIS provided. "
#                      "Define either METIS_DIR or METIS_LIB")
#endif(METIS_LIB_PATH)

########################################
# HWLOC

# Manage env
if(DEFINED ENV{HWLOC_LIB})
  set(HWLOC_LIB $ENV{HWLOC_LIB} CACHE PATH "Location of HWLOC library")
endif()
if(DEFINED ENV{HWLOC_INC})
  set(HWLOC_INC $ENV{HWLOC_INC} CACHE PATH "Location of HWLOC include")
endif()
if(DEFINED ENV{HWLOC_DIR})
  set(HWLOC_DIR $ENV{HWLOC_DIR} CACHE PATH "Location of HWLOC directory")
endif()

# Include
if(HWLOC_LIB AND HWLOC_INC)
  message(STATUS "HWLOC_INC used : " ${HWLOC_INC})
  message(STATUS "HWLOC_LIB used : " ${HWLOC_LIB})
  include_directories(${HWLOC_INC})
  set(HWLOC_LIB_PATH CACHE PATH ${HWLOC_LIB} "Location of HWLOC library")
else(HWLOC_LIB AND HWLOC_INC)
  if(HWLOC_DIR)
    include_directories(${HWLOC_DIR}/include)
    message(STATUS "HWLOC_DIR used: " ${HWLOC_DIR})
    set(HWLOC_LIB_PATH ${HWLOC_DIR}/hwloc/.libs CACHE PATH 
      "Location of HWLOC library")
  # else(HWLOC_DIR)
  #   message(ERROR " You should provide both HWLOC_LIB and HWLOC_INC, "
  #     "or HWLOC_DIR.")
  endif(HWLOC_DIR)
endif(HWLOC_LIB AND HWLOC_INC)

# Search library

unset(HWLOC_LIBRARIES)

if(HWLOC_LIB_PATH)

  find_library(HWLOC_LIBRARIES NAMES "hwloc" PATHS ${HWLOC_LIB_PATH})
  
else()

  find_library(HWLOC_LIBRARIES NAMES "hwloc")
  
endif(HWLOC_LIB_PATH)

if(HWLOC_LIBRARIES)
  set(LIBS ${LIBS} ${HWLOC_LIBRARIES})
else(HWLOC_LIBRARIES)
  message(FATAL_ERROR "HWLOC NOT FOUND")
endif(HWLOC_LIBRARIES)


#if(HWLOC_LIB)
#  set(HWLOC_LIB_PATH CACHE PATH ${HWLOC_LIB} "Location of HWLOC library")
#else(HWLOC_LIB)
#  if(DEFINED ENV{HWLOC_LIB})
#    set(HWLOC_LIB_PATH CACHE PATH $ENV{HWLOC_LIB} "Location of HWLOC library")
#  else(DEFINED ENV{HWLOC_LIB})
#    set(HWLOC_DIR $ENV{HWLOC_DIR} CACHE PATH "Location of HWLOC directory")
#    include_directories(${HWLOC_DIR}/include)
#    message(STATUS "HWLOC Directory: " ${HWLOC_DIR})
#    set(HWLOC_LIB_PATH ${HWLOC_DIR}/hwloc/.libs CACHE PATH 
#      "Location of HWLOC library")
#  endif(DEFINED ENV{HWLOC_LIB})
#endif(HWLOC_LIB)
#
#if(HWLOC_LIB_PATH)
#  set(LIB_TO_FIND "hwloc")
#
#  find_library(HWLOC_LIBRARIES NAMES ${LIB_TO_FIND} PATHS ${HWLOC_LIB_PATH})
#
#  if(HWLOC_LIBRARIES)
#    set(LIBS ${LIBS} ${HWLOC_LIBRARIES})
#   #message(STATUS "Found ${HWLOC_LIBRARIES}")
#   #set(CMAKE_REQUIRED_LIBRARIES ${HWLOC_LIBRARIES})
#   #CHECK_FUNCTION_EXISTS(${FUNCTION_TO_FIND} HWLOC_FOUND)
#   #unset(CMAKE_REQUIRED_LIBRARIES)
#   #message(STATUS "CHECK_FUNCTION_EXISTS on ${FUNCTION_TO_FIND}"
#   #                " returned : ${HWLOC_FOUND}")
#
#   #if(HWLOC_FOUND)
#   #  set(LIBS ${LIBS} ${HWLOC_LIBRARIES})
#   #else(HWLOC_FOUND)
#   #  message(ERROR "${FUNCTION_TO_FIND} can not be found in ${LIB_TO_FIND}")
#   #endif(HWLOC_FOUND)
#  else(HWLOC_LIBRARIES)
#    message(FATAL_ERROR "${HWLOC_LIB_PATH} does not contain lib${LIB_TO_FIND}")
#  endif(HWLOC_LIBRARIES)
#
#  unset(LIB_TO_FIND)
#else(HWLOC_LIB_PATH)
#  message(FATAL_ERROR "No Path for HWLOC provided. "
#                      "Define either HWLOC_DIR or HWLOC_LIB")
#endif(HWLOC_LIB_PATH)

########################################
# if (${BLAS_LAPACK} MATCHES "MKL")

#   # MKL library

#   if(${USE_KNL} MATCHES ON)

#     find_library(MKL_IFACE_LIBRARY mkl_intel_lp64 PATHS ${MKL_LIBS} ${MKL_DIR}/lib/intel64)
#     find_library(MKL_THREAD_LIBRARY mkl_intel_thread PATHS ${MKL_ROOT}/lib/intel64)
    
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
#     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")

#     target_link_libraries(spllt_test iomp5)

#   else()

#     set(MKL_LIBS $ENV{MKL_LIBS})
#     set(MKL_ROOT "/opt/intel/mkl" CACHE PATH "Location of MKL Library")
#     set(MKL_DIR $ENV{MKLROOT})

#     find_library(MKL_IFACE_LIBRARY mkl_gf_lp64 PATHS ${MKL_LIBS} ${MKL_ROOT}/lib/intel64 ${MKL_DIR}/lib/intel64)
#     # find_library(MKL_THREAD_LIBRARY mkl_gnu_thread PATHS ${MKL_ROOT}/lib/intel64)
#     find_library(MKL_SEQ_LIBRARY mkl_sequential PATHS ${MKL_LIBS} ${MKL_DEFAULT_ROOT}/lib/intel64 ${MKL_DIR}/lib/intel64)
#     find_library(MKL_CORE_LIBRARY mkl_core PATHS ${MKL_LIBS} ${MKL_ROOT}/lib/intel64 ${MKL_DIR}/lib/intel64) 
#     # set(MKL_LIBRARIES ${MKL_IFACE_LIBRARY} ${MKL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY})

#     MESSAGE( STATUS "MKL Root: " ${MKL_DIR} )
#     MESSAGE( STATUS "MKL MKL_IFACE_LIBRARY: " ${MKL_IFACE_LIBRARY} )

#   endif()

#   set(BLAS_LAPACK_LIBRARIES ${MKL_IFACE_LIBRARY} ${MKL_SEQ_LIBRARY} ${MKL_CORE_LIBRARY})
#   MESSAGE( STATUS "BLAS and LAPACK Libraries: " ${BLAS_LAPACK_LIBRARIES} )

# else()

#   # Reference BLAS libraries

#   find_library(BLAS_LIBRARY blas)
#   find_library(LAPACK_LIBRARY lapack)
#   set(BLAS_LAPACK_LIBRARIES ${BLAS_LIBRARY} ${LAPACK_LIBRARY})

# endif()

########################################
# BLAS
set(LBLAS "" CACHE STRING "BLAS library")

# If LBLAS not set, use environement variable
if(LBLAS)
  set(BLAS_LIBRARIES ${LBLAS})
elseif(DEFINED ENV{BLAS_LIB})
  set(BLAS_LIBRARIES $ENV{BLAS_LIB})
endif()

if(DEFINED BLAS_LIBRARIES)
  message(STATUS "BLAS : ${BLAS_LIBRARIES}")
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
  CHECK_FORTRAN_FUNCTION_EXISTS("dgemm" BLAS_FOUND)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(NOT BLAS_FOUND)
    message(ERROR "User supplied BLAS is NOT working")
  endif()
else()
  find_package(BLAS)
endif()

if(BLAS_FOUND)
  set(LIBS ${LIBS} ${BLAS_LIBRARIES})
else()
  message(FATAL_ERROR "BLAS NOT found")
endif (BLAS_FOUND)

########################################
# LAPACK
set(LLAPACK "" CACHE STRING "LAPACK library")
# If LBLAS not set, use environement variable 
if(LLAPACK)
  set(LAPACK_LIBRARIES ${LLAPACK})
elseif(DEFINED ENV{LAPACK_LIB})
  set(LAPACK_LIBRARIES $ENV{LAPACK_LIB})
endif()
if(DEFINED LAPACK_LIBRARIES)
  set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  CHECK_FORTRAN_FUNCTION_EXISTS("dpotrf" LAPACK_FOUND)
  unset(CMAKE_REQUIRED_LIBRARIES)
  if(NOT LAPACK_FOUND)
    message(ERROR "User supplied LAPACK is NOT working")
  endif()

else()
  find_package(LAPACK)
endif()

if (LAPACK_FOUND)
  set(LIBS ${LIBS} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
else()
  message(FATAL_ERROR "LAPACK library NOT found")
endif (LAPACK_FOUND)

########################################
# std C++ (needed for linking with ssids)
set(LIBS ${LIBS} stdc++)

# Add sources files
add_subdirectory(src)

################################################################################
# Test 

########################################
# Setup list of test files
set(test_files)

if(${RUNTIME} MATCHES "OMP")
  list(APPEND test_files "test/test_solve_phasis.F90")
endif()

foreach(test_file ${test_files})
  
  get_filename_component(test_file_name ${test_file} NAME_WE)
  add_executable(${test_file_name} ${test_file})

  ########################################
  # Link test drivers with external libraries
  target_include_directories(${test_file_name} PUBLIC ${SPRAL_DIR})
  target_include_directories(${test_file_name} PUBLIC ${CMAKE_BINARY_DIR}/src)
  target_link_libraries(${test_file_name} ${LIBS})

endforeach()


################################################################################
# Test drivers

########################################
# Setup list of test drivers
set(test_drivers)

#list(APPEND test_drivers "drivers/spllt_test.F90")

if(${RUNTIME} MATCHES "OMP")
  list(APPEND test_drivers "drivers/spllt_omp.F90")
  list(APPEND test_drivers "drivers/spllt_new_test.F90")
  list(APPEND test_drivers "drivers/spllt_omp_bench.F90")
endif()


option(C_IFACE_TEST "Compile the C examples" OFF)
if(C_IFACE_TEST MATCHES "ON")
  MESSAGE( STATUS "C Interface test: " ${C_IFACE_TEST} )

  include_directories("include/")
  list(APPEND c_iface_test "example/C/simple.c")

  if(DEFINED ENV{CPALAMEM_INC})
    set(CPALAMEM_LIB $ENV{CPALAMEM_LIB} CACHE PATH "Location of CPALAMEM library")
  endif()
  if(DEFINED ENV{CPALAMEM_INC})
    set(CPALAMEM_INC $ENV{CPALAMEM_INC} CACHE PATH "Location of CPALAMEM include")
  endif()

  if(CPALAMEM_INC AND CPALAMEM_LIB)
    message(STATUS "CPALAMEM_INC used : " ${CPALAMEM_INC})
    message(STATUS "CPALAMEM_LIB used : " ${CPALAMEM_LIB})
    include_directories(${CPALAMEM_INC})

    set(LIB_TO_FIND "cpalamem_core")

    find_library(CPALAMEM_LIBRARIES NAMES ${LIB_TO_FIND} PATHS ${CPALAMEM_LIB})

    if(CPALAMEM_LIBRARIES)
      set(LIBS ${LIBS} ${CPALAMEM_LIBRARIES})
    else(CPALAMEM_LIBRARIES)
      message(FATAL_ERROR "${CPALAMEM_LIB_PATH} does not contain lib${LIB_TO_FIND}")
    endif(CPALAMEM_LIBRARIES)

    unset(LIB_TO_FIND)
    list(APPEND c_iface_test "example/C/cpalamem_simple.c")
  else()
    message(WARNING "Can not set CPALAMEM.")
  endif()
endif()

option(PARDISO_TEST_DRIVER "Compile the Pardiso test driver" OFF)
if(PARDISO_TEST_DRIVER MATCHES "ON")

  # Manage env
  if(DEFINED ENV{PARDISO_INC})
    set(PARDISO_INC $ENV{PARDISO_INC} CACHE PATH "Location of PARDISO include")
  endif()

  # Include
  if(LPARDISO AND PARDISO_INC)
    message(STATUS "PARDISO_INC used  : " ${PARDISO_INC})
    message(STATUS "LPARDISO used     : " ${LPARDISO})
    include_directories(${PARDISO_INC})
    set(PARDISO_LIBRARIES ${LPARDISO})
  else(LPARDISO AND PARDISO_INC)
    message(ERROR " You should provide both LPARDISO and PARDISO_INC")
  endif(LPARDISO AND PARDISO_INC)

  # Search library
  if(PARDISO_LIBRARIES)
   #set(CMAKE_REQUIRED_LIBRARIES ${PARDISO_LIBRARIES})
   #message(STATUS "Search in " ${CMAKE_REQUIRED_LIBRARIES})
   #CHECK_FORTRAN_FUNCTION_EXISTS("pardiso" PARDISO_FOUND)
   #unset(CMAKE_REQUIRED_LIBRARIES)
   #if(NOT PARDISO_FOUND)
   #  message(ERROR " User supplied PARDISO library is NOT working")
   #else()
      set(LIBS ${LIBS} ${PARDISO_LIBRARIES})
      list(APPEND test_drivers "drivers/pardiso/pardiso_test.F90")
   #endif()
  else(PARDISO_LIBRARIES)
    message(FATAL_ERROR " No Path for PARDISO library provided. "
      "Define both LPARDISO and PARDISO_INC.")
  endif(PARDISO_LIBRARIES)

endif()

###################################
# Set the library installation
#set(PROJECT_VERSION "1.0.0")
#set_target_properties(spllt PROPERTIES VERSION ${PROJECT_VERSION})
#include(GNUInstallDirs)
#install(TARGETS spllt
#    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

foreach(test_driver ${test_drivers})
  
  get_filename_component(test_driver_name ${test_driver} NAME_WE)
  add_executable(${test_driver_name} ${test_driver})

  ########################################
  # Link test drivers with external libraries
  target_include_directories(${test_driver_name} PUBLIC ${SPRAL_DIR})
  target_include_directories(${test_driver_name} PUBLIC ${CMAKE_BINARY_DIR}/src)
  target_link_libraries(${test_driver_name} ${LIBS})

endforeach()

foreach(c_iface_test ${c_iface_test})
  
  get_filename_component(c_iface_test_name ${c_iface_test} NAME_WE)
  add_executable(${c_iface_test_name} ${c_iface_test})

  ########################################
  # Link test drivers with external libraries
  target_include_directories(${test_driver_name} PUBLIC ${SPRAL_DIR})
  target_include_directories(${test_driver_name} PUBLIC ${CMAKE_BINARY_DIR}/src)
 #target_include_directories(${test_driver_name} PUBLIC ${CMAKE_BINARY_DIR}/interfaces)
  target_link_libraries(${c_iface_test_name} ${LIBS})

endforeach()

################################################################################
# Print summary

message("-----------------------------------------")
message("Configuration of SpLLT done.")
message(" ")
message(
  "Build type: ${CMAKE_BUILD_TYPE}\n"
  "Fortran compiler: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})\n"
  "Fortran compiler flags: ${CMAKE_Fortran_FLAGS_ALL}\n"
  "C compiler:  ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})\n"
  "C compiler flags: ${CMAKE_Fortran_FLAGS_ALL}")
message("")
message("BLAS libraries: ${BLAS_LIBRARIES}")
message("LAPACK libraries: ${LAPACK_LIBRARIES}")
message("")
if(${RUNTIME} MATCHES "StarPU")
  message("Runtime system: StarPU")
  message("StarPU libraries: ${STARPU_LIBRARIES}")
elseif(${RUNTIME} MATCHES "Parsec")
  message("Runtime system: PaRSEC")
  message("Parsec libraries: ${PARSEC_LIBRARIES}")
elseif(${RUNTIME} MATCHES "OMP")
  message("Runtime system: OpenMP")
endif()
message("-----------------------------------------")
