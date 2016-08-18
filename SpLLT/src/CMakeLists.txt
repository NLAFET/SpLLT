# get source files
file(
  GLOB 
  sources 
  *.F90
  *.c)

if(${RUNTIME} MATCHES "StarPU")

  file(
    GLOB 
    starpu_sources 
    StarPU/*.F90
    StarPU/*.c)

  set(sources ${sources} ${starpu_sources})

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

set(HSLPACK_DIR $ENV{HSLPACKDIR})

set(HSL_SRC hsl_ma87d.F90 
  ${HSLPACK_DIR}/hsl_mc78/hsl_mc78i.f90 
  ${HSLPACK_DIR}/hsl_mc34/hsl_mc34d.f90
  ${HSLPACK_DIR}/hsl_zd11/hsl_zd11d.f90
  ${HSLPACK_DIR}/hsl_mc68/hsl_mc68i.f90
  ${HSLPACK_DIR}/hsl_zb01/hsl_zb01i.f90
  ${HSLPACK_DIR}/hsl_mc34/hsl_mc34d.f90
  ${HSLPACK_DIR}/hsl_mc56/hsl_mc56d.f90
  ${HSLPACK_DIR}/hsl_fa14/hsl_fa14d.f90
  ${HSLPACK_DIR}/archive/hsl_zd13d.f90
  ${HSLPACK_DIR}/hsl_mc69/hsl_mc69d.f90
  ${HSLPACK_DIR}/fa14/fa14d.f
  ${HSLPACK_DIR}/mc34/mc34d.f
  ${HSLPACK_DIR}/mc47/mc47d.f
  ${HSLPACK_DIR}/ym11/ym11d.f
  ${HSLPACK_DIR}/mc59/mc59d.f
  ${HSLPACK_DIR}/mc54/mc54d.f
  ${HSLPACK_DIR}/mc56/mc56d.f
  ${HSLPACK_DIR}/mc56/mc56i.f
  ${HSLPACK_DIR}/fa14/fa14d.f
  )

# MESSAGE( STATUS "source files:         " ${HSLSRC} )

# add_library(spllt ${sources} ${HSL_SRC})
# MESSAGE( STATUS "CUDA found:         " ${CUDA_FOUND} )
if(${USE_GPU} MATCHES ON)
  cuda_add_library(spllt_cuda ${cuda_sources})
endif()

add_library(spllt ${sources} ${HSL_SRC})
