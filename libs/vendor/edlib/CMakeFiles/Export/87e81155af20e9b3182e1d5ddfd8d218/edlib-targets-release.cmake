#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "edlib::edlib" for configuration "Release"
set_property(TARGET edlib::edlib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(edlib::edlib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libedlib.a"
  )

list(APPEND _cmake_import_check_targets edlib::edlib )
list(APPEND _cmake_import_check_files_for_edlib::edlib "${_IMPORT_PREFIX}/lib/libedlib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
