# get source files
file(GLOB sources      *.F90)

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
  ${HSLPACK_DIR}/fa14/fa14d.f
  ${HSLPACK_DIR}/mc34/mc34d.f
  ${HSLPACK_DIR}/mc47/mc47d.f
  ${HSLPACK_DIR}/ym11/ym11d.f
  ${HSLPACK_DIR}/mc59/mc59d.f
  ${HSLPACK_DIR}/mc54/mc54d.f
  )

# MESSAGE( STATUS "source files:         " ${HSLSRC} )

add_library(spllt ${sources} ${HSL_SRC})
