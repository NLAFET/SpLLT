# get source files

file(
  GLOB 
  sources 
  *.F90
  *.c)

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


# MESSAGE( STATUS "CUDA found:         " ${CUDA_FOUND} )
if(${USE_GPU} MATCHES ON)
  cuda_add_library(spllt_cuda ${cuda_sources})
endif()

if(WITH_TRACE MATCHES ON)
  add_library(spllt_trace ${sources} ${sources_iface})
else()
  add_library(spllt ${sources} ${sources_iface})
endif()
