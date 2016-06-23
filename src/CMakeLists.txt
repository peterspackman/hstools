set(API_SRC
    hs.F90
    api.F90
    h5file.F90
)
add_library(tontoapi SHARED ${API_SRC})
target_link_libraries(tontoapi tonto ${HDF5_LIBRARIES} ${FGSL_LIBRARIES} ${LAPACK_LIBRARIES})
