# DO NOT EDIT THIS FILE! THIS FILE IS AUTOMATICALLY GENERATED BY xml2cmake
PROJECT(STAR)

INCLUDE(${VISIT_SOURCE_DIR}/CMake/PluginMacros.cmake)

SET(COMMON_SOURCES
STARPluginInfo.C
STARCommonPluginInfo.C
)

SET(LIBI_SOURCES 
STARPluginInfo.C
)

SET(LIBM_SOURCES
STARMDServerPluginInfo.C
${COMMON_SOURCES}
avtSTARFileFormat.C
ConfigFileReader.C
MultiresFileReader.C
MultiresGrid.C
MultiresMultifileReader.C
ResolutionMap.C
)

SET(LIBE_SOURCES
STAREnginePluginInfo.C
${COMMON_SOURCES}
avtSTARFileFormat.C
ConfigFileReader.C
MultiresFileReader.C
MultiresGrid.C
MultiresMultifileReader.C
ResolutionMap.C
)

INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}
${VISIT_COMMON_INCLUDES}
${VISIT_INCLUDE_DIR}/avt/DBAtts/MetaData
${VISIT_INCLUDE_DIR}/avt/DBAtts/SIL
${VISIT_INCLUDE_DIR}/avt/Database/Database
${VISIT_INCLUDE_DIR}/avt/Database/Formats
${VISIT_INCLUDE_DIR}/avt/Database/Ghost
${VISIT_INCLUDE_DIR}/avt/FileWriter
${VISIT_INCLUDE_DIR}/avt/Filters
${VISIT_INCLUDE_DIR}/avt/MIR/Base
${VISIT_INCLUDE_DIR}/avt/MIR/Tet
${VISIT_INCLUDE_DIR}/avt/MIR/Zoo
${VISIT_INCLUDE_DIR}/avt/Math
${VISIT_INCLUDE_DIR}/avt/Pipeline/AbstractFilters
${VISIT_INCLUDE_DIR}/avt/Pipeline/Data
${VISIT_INCLUDE_DIR}/avt/Pipeline/Pipeline
${VISIT_INCLUDE_DIR}/avt/Pipeline/Sinks
${VISIT_INCLUDE_DIR}/avt/Pipeline/Sources
${VISIT_INCLUDE_DIR}/avt/VisWindow/VisWindow
${VISIT_INCLUDE_DIR}/visit_vtk/full
${VISIT_INCLUDE_DIR}/visit_vtk/lightweight
${EAVL_INCLUDE_DIR} 
${VTK_INCLUDE_DIRS} 
)

ADD_DEFINITIONS(-DSTAR_VISIT)

LINK_DIRECTORIES(${VISIT_LIBRARY_DIR} ${EAVL_LIBRARY_DIR} ${VTK_LIBRARY_DIRS} )

ADD_LIBRARY(ISTARDatabase ${LIBI_SOURCES})
TARGET_LINK_LIBRARIES(ISTARDatabase visitcommon)
SET(INSTALLTARGETS ISTARDatabase)

IF(NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)
    ADD_LIBRARY(MSTARDatabase ${LIBM_SOURCES}    )
    TARGET_LINK_LIBRARIES(MSTARDatabase visitcommon avtdbatts avtdatabase_ser )
    ADD_TARGET_DEFINITIONS(MSTARDatabase MDSERVER)
    SET(INSTALLTARGETS ${INSTALLTARGETS} MSTARDatabase)
ENDIF(NOT VISIT_ENGINE_ONLY AND NOT VISIT_DBIO_ONLY)

ADD_LIBRARY(ESTARDatabase_ser ${LIBE_SOURCES})
TARGET_LINK_LIBRARIES(ESTARDatabase_ser visitcommon avtdatabase_ser avtpipeline_ser )
ADD_TARGET_DEFINITIONS(ESTARDatabase_ser ENGINE)
SET(INSTALLTARGETS ${INSTALLTARGETS} ESTARDatabase_ser)

IF(VISIT_PARALLEL)
    ADD_PARALLEL_LIBRARY(ESTARDatabase_par ${LIBE_SOURCES})
    TARGET_LINK_LIBRARIES(ESTARDatabase_par visitcommon avtdatabase_par avtpipeline_par )
    ADD_TARGET_DEFINITIONS(ESTARDatabase_par ENGINE)
    SET(INSTALLTARGETS ${INSTALLTARGETS} ESTARDatabase_par)
ENDIF(VISIT_PARALLEL)

VISIT_INSTALL_DATABASE_PLUGINS(${INSTALLTARGETS})
VISIT_PLUGIN_TARGET_RTOD(databases ${INSTALLTARGETS})
VISIT_PLUGIN_TARGET_FOLDER(databases STAR ${INSTALLTARGETS})

