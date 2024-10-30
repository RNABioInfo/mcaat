#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "edlib::edlib" for configuration "Debug"
set_property(TARGET edlib::edlib APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(edlib::edlib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libedlib.a"
  )

list(APPEND _cmake_import_check_targets edlib::edlib )
list(APPEND _cmake_import_check_files_for_edlib::edlib "${_IMPORT_PREFIX}/lib/libedlib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
