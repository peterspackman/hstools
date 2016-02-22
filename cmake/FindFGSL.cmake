# - Try to find FGSL
# Once done this will define
#  FGSL_FOUND - System has LibXml2
#  FGSL_INCLUDE_DIRS - The LibXml2 include directories
#  FGSL_LIBRARIES - The libraries needed to use LibXml2
#  FGSL_DEFINITIONS - Compiler switches required for using LibXml2

find_package(PkgConfig)
pkg_check_modules(PC_FGSL fgsl)
set(FGSL_DEFINITIONS ${PC_FGSL_CFLAGS_OTHER})

find_path(FGSL_INCLUDE_DIR fgsl/fgsl.mod
  HINTS ${PC_FGSL_INCLUDEDIR} ${PC_FGSL_INCLUDE_DIRS}
          PATH_SUFFIXES fgsl )

find_library(FGSL_LIBRARY NAMES fgsl
   HINTS ${PC_FGSL_LIBDIR} ${PC_FGSL_LIBRARY_DIRS} )

set(FGSL_LIBRARIES ${FGSL_LIBRARY} )
set(FGSL_INCLUDE_DIRS ${FGSL_INCLUDE_DIR}/fgsl )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FGSL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(FGSL  DEFAULT_MSG
                                  FGSL_LIBRARY FGSL_INCLUDE_DIR)

mark_as_advanced(FGSL_INCLUDE_DIR FGSL_LIBRARY )
