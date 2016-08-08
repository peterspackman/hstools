set(API_SRC
    hs.F90
    api.F90
    h5file.F90
)
add_library(tontoapi SHARED ${API_SRC})
target_link_libraries(tontoapi tonto ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES} ${LAPACK_LIBRARIES})
