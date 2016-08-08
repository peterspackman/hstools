##################################################
# Define the project and associated dependencies
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
PROJECT(hstools LANGUAGES Fortran C VERSION 0.10.0)

# Set the hs version

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})
set(CMAKE_INSTALL_PREFIX CACHE PATH "/usr/local/")
set(HDF5_USER_STATIC_LIBRARIES TRUE)

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
add_subdirectory(${CMAKE_SOURCE_DIR}/tonto
    ${TONTO_BUILD_DIR}
    EXCLUDE_FROM_ALL)

find_package(LAPACK REQUIRED)

################ HDF5 library #################
find_package(HDF5 COMPONENTS Fortran_HL Fortran C HL)

include_directories("${HDF5_INCLUDE_DIRS}")
include_directories("${TONTO_BUILD_DIR}")
add_executable(hs ${SRC})
target_link_libraries(hs tonto ${HDF5_Fortran_LIBRARIES}
    ${HDF5_Fortran_HL_LIBRARIES} ${LAPACK_LIBRARIES})
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src")
install(TARGETS hs DESTINATION bin)

enable_testing()
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/tests")

find_package(PythonInterp 3)

if(PythonInterp_FOUND) 
    cmake_policy(SET CMP0026 OLD)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(DEPS "${CMAKE_SOURCE_DIR}/hstools/__init__.py" hs)
    set(BUILD_TIMESTAMP "${CMAKE_CURRENT_BINARY_DIR}/timestamp")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})
    add_custom_command(OUTPUT ${BUILD_TIMESTAMP}
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build
        COMMAND ${CMAKE_COMMAND} -E touch ${BUILD_TIMESTAMP}
        DEPENDS ${DEPS})

    add_custom_target(py_install DEPENDS ${BUILD_TIMESTAMP})
    install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install)"
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
else()
    message(STATUS "Python not found")
endif()
