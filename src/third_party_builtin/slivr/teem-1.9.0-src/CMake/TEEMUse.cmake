IF(NOT TEEM_FOUND)
  MESSAGE(FATAL_ERROR "Something went wrong. You are including TEEMUse.cmake but TEEM was not found")
ENDIF(NOT TEEM_FOUND)

# Make TEEM easier to use
INCLUDE_DIRECTORIES(${TEEM_INCLUDE_DIRS})
LINK_DIRECTORIES(${TEEM_LIBRARY_DIRS})

# Load the compiler settings used for TEEM.
IF(TEEM_BUILD_SETTINGS_FILE)
  INCLUDE(CMakeImportBuildSettings)
  CMAKE_IMPORT_BUILD_SETTINGS(${TEEM_BUILD_SETTINGS_FILE})
ENDIF(TEEM_BUILD_SETTINGS_FILE)

# Add compiler flags needed to use TEEM.
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${TEEM_REQUIRED_C_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${TEEM_REQUIRED_EXE_LINKER_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${TEEM_REQUIRED_SHARED_LINKER_FLAGS}")
SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${TEEM_REQUIRED_MODULE_LINKER_FLAGS}")

