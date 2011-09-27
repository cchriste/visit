include(${VISIT_SOURCE_DIR}/config-site/windows.cmake)

# can provide overrides for windows.cmake settings below

# disable some warnings
# 4244 conversion double to float etc
# 4305 truncation from 'double' to 'float', etc
# 4800 'int' forcing value to bool (performance warning)
ADD_DEFINITIONS(/wd4244 /wd4305 /wd4800)

VISIT_OPTION_DEFAULT(VISIT_INSTALL_THIRD_PARTY ON TYPE BOOL)
VISIT_OPTION_DEFAULT(VISIT_PARALLEL ON TYPE BOOL)

#VISIT_OPTION_DEFAULT(VISIT_MPI_LIBS msmpi TYPE STRING)
##Enclose the include and lib paths in escape quotes in case of spaces.
#VISIT_OPTION_DEFAULT(VISIT_MPI_CXX_FLAGS "/I \"$ENV{CCP_INC}\"" TYPE STRING)
#VISIT_OPTION_DEFAULT(VISIT_MPI_C_FLAGS   "/I \"$ENV{CCP_INC}\"" TYPE STRING)
#IF(CMAKE_CL_64)
#    VISIT_OPTION_DEFAULT(VISIT_MPI_LD_FLAGS  "/LIBPATH:\"$ENV{CCP_LIB64}\"" TYPE STRING)
#ELSE(CMAKE_CL_64)
#    VISIT_OPTION_DEFAULT(VISIT_MPI_LD_FLAGS  "/LIBPATH:\"$ENV{CCP_LIB32}\"" TYPE STRING)
#ENDIF(CMAKE_CL_64)

