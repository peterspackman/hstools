set(API_SRC
    hs.F90
    api.F90
    h5file.F90
)
add_library(hsapi SHARED ${API_SRC})
target_link_libraries(hsapi tonto ${HDF5_LIBRARIES} ${FGSL_LIBRARIES} ${LAPACK_LIBRARIES})
