#*****************************************************************************
#
# Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
# Produced at the Lawrence Livermore National Laboratory
# LLNL-CODE-442911
# All rights reserved.
#
# This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
# full copyright notice is contained in the file COPYRIGHT located at the root
# of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
#
# Redistribution  and  use  in  source  and  binary  forms,  with  or  without
# modification, are permitted provided that the following conditions are met:
#
#  - Redistributions of  source code must  retain the above  copyright notice,
#    this list of conditions and the disclaimer below.
#  - Redistributions in binary form must reproduce the above copyright notice,
#    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
#    documentation and/or other materials provided with the distribution.
#  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
#    be used to endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
# ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
# LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
# DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
# SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
# CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
# LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
# OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#
# Modifications:
#   Kathleen Biagas, Thu June 14 15:48:01 MST 2012
#   Add support for pyside 1.1.1 on windows, which no longer has a separate 
#   generatorrunner module.  For PYSIDE_ADD_MODULE function, change the 
#   target's output location on windows to the exe dir, explanation for this
#   is within the code below.  Change gen_pathsep on windows to "\;".
#
#   Gunther H. Weber, Thu Nov 21 18:28:56 PST 2013
#   Allow user to set Shiboken_FOUND and PySide_FOUND to set up relevant
#   CMake variables manually (necessary for MacPorts PySide install).
#   Use ${PYSIDE_PYTHONPATH} instead of ${VISIT_PYSIDE_DIR}/lib/site-packages
#   to support using system PySide.
#
#   Kathleen Biagas, Wed Dec 18 08:05:38 PST 2013
#   Exclude '.svn' from INSTALL command, probably a rare case, but a
#   developer did run into an issue.
#
#****************************************************************************/

INCLUDE(${VISIT_SOURCE_DIR}/CMake/SetUpThirdParty.cmake)

IF(VISIT_PYSIDE_DIR)

    #  Find PySide 
    set(CMAKE_PREFIX_PATH ${VISIT_PYSIDE_DIR}/lib/cmake/ ${CMAKE_PREFIX_PATH})
    set(CMAKE_LIBRARY_PATH ${VISIT_PYSIDE_DIR}/lib ${CMAKE_LIBRARY_PATH})

    if(NOT Shiboken_FOUND)
        find_package(Shiboken 1.1.1)
    endif(NOT Shiboken_FOUND)
    if (NOT PySide_FOUND)
        find_package(PySide 1.1.1)
    endif(NOT PySide_FOUND)
    IF(Shiboken_FOUND)
        SET(GENERATORRUNNER_BINARY ${SHIBOKEN_BINARY})
    ENDIF(Shiboken_FOUND)


ENDIF(VISIT_PYSIDE_DIR)

IF(NOT PySide_FOUND OR NOT Shiboken_FOUND)
    #If we dont have shiboken, force pyside off
    MESSAGE(STATUS "PySide NOT found")
    SET(PySide_FOUND 0)
ENDIF (NOT PySide_FOUND OR NOT Shiboken_FOUND)


IF(PySide_FOUND)
    SET_UP_THIRD_PARTY(PYSIDE lib include
          pyside-python${PYTHON_VERSION} shiboken-python${PYTHON_VERSION})
    # The PySide module is symlinked into the python install VisIt uses for 
    # dev builds.  For 'make install' and 'make package' we need to actually 
    # install the PySide SOs.
    SET(PYSIDE_MODULE_SRC  ${PYSIDE_PYTHONPATH}/PySide/)
    IF(UNIX)
        #SET(PYSIDE_MODULE_INSTALLED_DIR ${VISIT_INSTALLED_VERSION_LIB}/python/lib/python${PYTHON_VERSION}/site-packages/PySide/)
        SET(PYSIDE_MODULE_INSTALLED_DIR ${VISIT_INSTALLED_VERSION_LIB}/site-packages/PySide/)
    ELSE(UNIX)
        #SET(PYSIDE_MODULE_INSTALLED_DIR ${VISIT_INSTALLED_VERSION_LIB}/python/Lib/site-packages/PySide/)
        SET(PYSIDE_MODULE_INSTALLED_DIR ${VISIT_INSTALLED_VERSION_LIB}/site-packages/PySide/)
    ENDIF(UNIX)

    FILE(GLOB pysidelibs ${PYSIDE_MODULE_SRC}/*)
    INSTALL(FILES ${pysidelibs}
            DESTINATION ${PYSIDE_MODULE_INSTALLED_DIR}
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
                        GROUP_READ GROUP_WRITE GROUP_EXECUTE 
                        WORLD_READ WORLD_EXECUTE
            CONFIGURATIONS "" None Debug Release RelWithDebInfo MinSizeRel
            PATTERN ".svn" EXCLUDE
            )
ENDIF(PySide_FOUND)




#****************************************************************************
# PYSIDE_ADD_MODULE
# Defines a new PySide module and creates a dependent generator target and
# distutils setup call.
#****************************************************************************
FUNCTION(PYSIDE_ADD_MODULE module_name
                           dest_dir
                           mod_sources
                           mod_gen_sources
                           mod_gen_include_paths
                           mod_gen_link_libs
                           mod_gen_global
                           mod_gen_typesystem)

MESSAGE(STATUS "Configuring PySide module: ${module_name}")


PYSIDE_ADD_GENERATOR_TARGET("${module_name}_gen"
                            ${mod_gen_sources}
                            ${mod_gen_include_paths}
                            ${mod_gen_global}
                            ${mod_gen_typesystem})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${SHIBOKEN_INCLUDE_DIR}
                    ${SHIBOKEN_PYTHON_INCLUDE_DIR}
                    ${PYSIDE_INCLUDE_DIR}
                    ${CMAKE_CURRENT_BINARY_DIR}/${module_name}
                    ${${mod_gen_include_paths}})

add_library(${module_name} MODULE ${${mod_sources}} ${${mod_gen_sources}})

SET_TARGET_PROPERTIES(${module_name} PROPERTIES PREFIX "")

IF(UNIX)
    SET_TARGET_PROPERTIES(${module_name} PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${VISIT_LIBRARY_DIR}/${dest_dir})
ELSE(UNIX)
    # Use a different LIBARARY_OUTPUT_DIRECTORY on windows so visit -cli -pyside
    # will work from the build directory.  Here's why:
    #
    # If VISIT_LIBRARY_DIR/dest_dir is used, the files are output to
    # lib/Release/visit/site-packages/Release.  This causes problems
    # with running visit -cli -pyside from the build directory, because
    # the extra 'Release' at the end of the path is not part of sys.path.
    # Copying the output files to the directory one below solves the problem
    # of running from the build, but the INSTALL target fails because 
    # CMake only translates one of the $(Configuration) macros in the path:
    # "lib/$(Configuration)/visit/site-packages/$(Configuration)", so the file
    # "lib/$(Configuration)/visit/site-packages/Release" cannot be
    # found and the INSTALL target fails.
    # 
    # using VISIT_BINARY_DIR/exe works for running from the build dir because
    # the exe path is part of sys.path and the .pyd files are found and loaded
    # correctly.  It also works for the INSTALL target, because the files
    # are copied from the exe dir to the install-dir/lib/site-packages/visit.
    # 
    SET_TARGET_PROPERTIES(${module_name} PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${VISIT_BINARY_DIR}/exe)
    SET_TARGET_PROPERTIES(${module_name} PROPERTIES SUFFIX ".pyd")
ENDIF(UNIX)

target_link_libraries(${module_name}
                      ${SHIBOKEN_PYTHON_LIBRARIES}
                      ${SHIBOKEN_LIBRARY}
                      ${PYSIDE_LIBRARY}
                      ${${mod_gen_link_libs}})

add_dependencies(${module_name} "${module_name}_gen")

VISIT_INSTALL_TARGETS_RELATIVE(${dest_dir} ${module_name})

ENDFUNCTION(PYSIDE_ADD_MODULE)

#****************************************************************************
# PYSIDE_ADD_HYBRID_MODULE
# Defines a new PySide module and creates a dependent generator target and
# distutils setup call.
#****************************************************************************
FUNCTION(PYSIDE_ADD_HYBRID_MODULE module_name
                                  dest_dir
                                  mod_py_setup
                                  mod_py_sources
                                  mod_sources
                                  mod_gen_sources
                                  mod_gen_include_paths
                                  mod_gen_link_libs
                                  mod_gen_global
                                  mod_gen_typesystem)

MESSAGE(STATUS "Configuring PySide module: ${module_name}")

PYTHON_ADD_DISTUTILS_SETUP("${module_name}_py_setup"
                            ${dest_dir}
                            ${mod_py_setup}
                            ${mod_py_sources})

PYSIDE_ADD_MODULE(${module_name}
                  ${dest_dir}/${module_name}
                  ${mod_sources}
                  ${mod_gen_sources}
                  ${mod_gen_include_paths}
                  ${mod_gen_link_libs}
                  ${mod_gen_global}
                  ${mod_gen_typesystem})

add_dependencies(${module_name} "${module_name}_py_setup")
add_dependencies(${module_name} "${module_name}_gen")

ENDFUNCTION(PYSIDE_ADD_HYBRID_MODULE)


#****************************************************************************
# PYSIDE_ADD_GENERATOR_TARGET
# Defines a new PySide generator target.
#****************************************************************************
FUNCTION(PYSIDE_ADD_GENERATOR_TARGET 
            target_name gen_sources gen_include_paths gen_global gen_typesystem)

set(PYSIDE_GENERATOR_EXTRA_FLAGS --silent 
                                 --generator-set=shiboken 
                                 --enable-parent-ctor-heuristic 
                                 --enable-return-value-heuristic 
                                 --avoid-protected-hack)

#
# Create string for --include-paths arg from gen_include_paths list.
#

IF(WIN32)
    set(gen_pathsep "\;")
ELSE(WIN32)
    set(gen_pathsep ":")
ENDIF(WIN32)

SET(gen_include_paths_arg "${CMAKE_CURRENT_SOURCE_DIR}")

FOREACH(itm ${${gen_include_paths}})
    SET(gen_include_paths_arg "${gen_include_paths_arg}${gen_pathsep}${itm}")
ENDFOREACH(itm ${gen_include_paths})

add_custom_command(OUTPUT ${${gen_sources}}
                   COMMAND ${GENERATORRUNNER_BINARY} ${PYSIDE_GENERATOR_EXTRA_FLAGS}
                   ${gen_global}
                   --include-paths=${gen_include_paths_arg}
                   --typesystem-paths=${PYSIDE_TYPESYSTEMS}
                   --output-directory=${CMAKE_CURRENT_BINARY_DIR}
                   ${gen_typesystem}
                   DEPENDS  ${gen_global} ${gen_typesystem}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   COMMENT "Running generator for ${${gen_sources}}...")

add_custom_target(${target_name} DEPENDS ${${gen_sources}})
ENDFUNCTION(PYSIDE_ADD_GENERATOR_TARGET)

