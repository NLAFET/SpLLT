# get source files

file(
  GLOB 
  sources 
  *.F90
  *.c)

file(
  GLOB 
  to_remove
  hsl_ma87d.F90)

list(REMOVE_ITEM sources ${to_remove})

file(
  GLOB 
  sources_iface
  ../interfaces/C/*.F90)

if(${RUNTIME} MATCHES "StarPU")

  file(
    GLOB 
    starpu_sources 
    StarPU/*.F90
    StarPU/*.c)

  set(sources ${sources} ${starpu_sources})

elseif(${RUNTIME} MATCHES "Parsec")

  # adds headers to the include directories
  include_directories(${CMAKE_SOURCE_DIR}/src/PaRSEC)

  # include(${PARSEC_SRC_DIR}/cmake_modules/RulesJDF.cmake)
  include(${CMAKE_SOURCE_DIR}/cmake/Modules/RulesJDF.cmake)

  jdf_rules(generated_jdf "${CMAKE_SOURCE_DIR}/src/PaRSEC/factorize.jdf")
  jdf_rules(generated_jdf "${CMAKE_SOURCE_DIR}/src/PaRSEC/gather.jdf")

  # message(STATUS "Generated JDFs: ${generated_jdf}")

  file(
    GLOB 
    parsec_sources
    PaRSEC/*.F90
    PaRSEC/*.c)

  file(
    GLOB 
    to_remove
    PaRSEC/spllt_parsec_data.c)
  
  list(REMOVE_ITEM parsec_sources ${to_remove})

  set(sources ${sources} ${generated_jdf} ${parsec_sources})

  # message(STATUS "Parsec sources: ${sources}")

endif()

if(${USE_GPU} MATCHES ON)

  file(
    GLOB 
    cuda_sources 
    StarPU/*.cu)

  # set(sources ${sources} ${cuda_sources})
  # set(sources ${cuda_sources})
  
  # MESSAGE( STATUS "CUDA source files:         " ${sources} )
endif()

# add_library(hsl_ma87 hsl_ma87d.F90
#   $ENV{HSLPACKDIR}/hsl_mc78/hsl_mc78i.f90
#   $ENV{HSLPACKDIR}/hsl_mc34/hsl_mc34d.f90)

# if(MA87_TEST_DRIVER MATCHES ON)

set(HSLPACK_DIR $ENV{HSLPACK_DIR})

# set(HSL_SRC hsl_ma87d.F90
#   ${HSLPACK_DIR}/hsl_mc78/hsl_mc78i.f90
#   ${HSLPACK_DIR}/hsl_mc34/hsl_mc34d.f90
#   ${HSLPACK_DIR}/hsl_zd11/hsl_zd11d.f90
#   ${HSLPACK_DIR}/hsl_mc68/hsl_mc68i.f90
#   ${HSLPACK_DIR}/hsl_zb01/hsl_zb01i.f90
#   ${HSLPACK_DIR}/hsl_mc34/hsl_mc34d.f90
#   ${HSLPACK_DIR}/hsl_mc56/hsl_mc56d.f90
#   ${HSLPACK_DIR}/hsl_fa14/hsl_fa14d.f90
#   ${HSLPACK_DIR}/archive/hsl_zd13d.f90
#   ${HSLPACK_DIR}/hsl_mc69/hsl_mc69d.f90
#   ${HSLPACK_DIR}/fa14/fa14d.f
#   ${HSLPACK_DIR}/mc34/mc34d.f
#   ${HSLPACK_DIR}/mc47/mc47d.f
#   ${HSLPACK_DIR}/ym11/ym11d.f
#   ${HSLPACK_DIR}/mc59/mc59d.f
#   ${HSLPACK_DIR}/mc54/mc54d.f
#   ${HSLPACK_DIR}/mc56/mc56d.f
#   ${HSLPACK_DIR}/mc56/mc56i.f
#   ${HSLPACK_DIR}/fa14/fa14d.f
#   )

if(MA87_TEST_DRIVER MATCHES ON)
set(HSL_SRC hsl_ma87d.F90
  ${HSLPACK_DIR}/hsl_mc78/hsl_mc78i.f90
  ${HSLPACK_DIR}/hsl_mc34/hsl_mc34d.f90)
endif(MA87_TEST_DRIVER MATCHES ON)

# endif()
# MESSAGE( STATUS "source files:         " ${HSLSRC} )

# add_library(spllt ${sources} ${HSL_SRC})
# MESSAGE( STATUS "CUDA found:         " ${CUDA_FOUND} )
if(${USE_GPU} MATCHES ON)
  cuda_add_library(spllt_cuda ${cuda_sources})
endif()

if(MA87_TEST_DRIVER MATCHES ON)
  # MESSAGE( STATUS "Sources: " ${sources} )
  # MESSAGE( STATUS "HSL Sources: " ${HSL_SRC} )
 #add_library(spllt_lib ${sources} ${HSL_SRC} ${sources_iface})
  if(WITH_TRACE MATCHES ON)
    add_library(spllt_trace ${sources} ${HSL_SRC} ${sources_iface})
  else()
    add_library(spllt ${sources} ${HSL_SRC} ${sources_iface})
  endif()
else()
  if(WITH_TRACE MATCHES ON)
    add_library(spllt_trace ${sources} ${sources_iface})
  else()
    add_library(spllt ${sources} ${sources_iface})
  endif()
endif()
