##################################################
# Define the project and associated dependencies
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.1)
PROJECT(hs Fortran)

# Set the hs version
SET(VERSION 0.10.0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})
set(CMAKE_INSTALL_PREFIX CACHE PATH "/usr/local/")

# Set up defaults for this project
include(ExternalProject)

set(SRC
    src/hs.f90
    src/h5file.f90
    src/main.f90
)

enable_testing()

###############################################
################ Tonto library ################
###############################################
set(TONTO_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/tonto)
add_subdirectory(${CMAKE_SOURCE_DIR}/externals/tonto
    ${TONTO_BUILD_DIR}
    EXCLUDE_FROM_ALL)

find_package(LAPACK REQUIRED)

###############################################
################ FGSL library #################
###############################################
find_package(FGSL REQUIRED)

if(NOT ${FGSL_FOUND})
    message("FGSL not found; will download and compile instead")
    SET(FGSL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fgsl-1.0.0)
    ExternalProject_Add(project_fgsl
        URL http://www.lrz.de/services/software/mathematik/gsl/fortran/fgsl-1.0.0.tar.gz
        PREFIX ${FGSL_PREFIX}
        CONFIGURE_COMMAND "<SOURCE_DIR>/configure"
        BUILD_COMMAND make
        INSTALL_COMMAND make install)

    ExternalProject_Get_Property(project_fgsl install_dir)
    add_dependencies(fgsl project_fgsl)
    set_property(TARGET fgsl PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libfgsl.a)
    add_library(fgsl STATIC IMPORTED)
    set(FGSL_INCLUDE_DIRS ${install_dir}/include)
    set(FGSL_LIBRARIES ${install_dir}/lib/libfgsl.a})
endif ()

###############################################
################ HDF5 library #################
###############################################
find_package(HDF5 COMPONENTS Fortran_HL Fortran C HL)

if(NOT ${HDF5_FOUND})
    message("HDF5 not found; will download and compile instead")
    SET(HDF5_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/hdf5)
    SET(HDF5_CONFIG "--enable-fortran --enable-fortran2003 --enable-hl --enable-static")
    ExternalProject_Add(project_hdf5
        URL http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.15-patch1.tar.gz
        PREFIX ${HDF5_PREFIX}
        CONFIGURE_COMMAND "<SOURCE_DIR>/configure ${HDF5_CONFIG}"
        BUILD_COMMAND make
        INSTALL_COMMAND make install)

    ExternalProject_Get_Property(project_hdf5 install_dir)
    add_dependencies(hdf5 project_hdf5)
    set_property(TARGET hdf5 PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libfgsl.a)
    add_library(hdf5 STATIC IMPORTED)
    set(HDF5_INCLUDE_DIRS ${install_dir}/include)
    set(HDF5_LIBRARIES "${install_dir}/lib/libhdf5.a}"
        "${install_dir}/lib/libhdf5_fortran.a}"
        "${install_dir}/lib/libhdf5_hl.a}"
        "${install_dir}/lib/libhdf5hl_fortran.a}"
    )

endif()

include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${FGSL_INCLUDE_DIRS}")
include_directories("${TONTO_BUILD_DIR}")
add_executable(hs ${SRC})
target_link_libraries(hs tonto ${HDF5_LIBRARIES} ${FGSL_LIBRARIES} ${LAPACK_LIBRARIES})
install(TARGETS hs DESTINATION bin)
