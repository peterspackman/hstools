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
    src/hs.F90
    src/h5file.F90
    src/main.f90
    src/api.F90
)


################ Tonto library ################
set(TONTO_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/tonto)
add_subdirectory(${CMAKE_SOURCE_DIR}/externals/tonto
    ${TONTO_BUILD_DIR}
    EXCLUDE_FROM_ALL)

find_package(LAPACK REQUIRED)

################ HDF5 library #################
find_package(HDF5 COMPONENTS Fortran_HL Fortran C HL)

include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${TONTO_BUILD_DIR}")
add_executable(hs ${SRC})
target_link_libraries(hs tonto ${HDF5_LIBRARIES} ${FGSL_LIBRARIES} ${LAPACK_LIBRARIES})
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src")
install(TARGETS hs DESTINATION bin)

enable_testing()
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/tests")
