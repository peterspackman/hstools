add_executable(h5test 
    ${CMAKE_CURRENT_SOURCE_DIR}/test_h5file.f90
    ${PROJECT_SOURCE_DIR}/src/h5file.F90)
target_link_libraries(h5test ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
add_test(NAME "H5File" COMMAND h5test)
