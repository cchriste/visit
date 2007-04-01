# Microsoft Developer Studio Project File - Name="PDBM" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=PDBM - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "PDBM.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "PDBM.mak" CFG="PDBM - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "PDBM - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "PDBM - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "PDBM - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "SILOM_EXPORTS" /YX /FD /c
# ADD CPP /nologo /G6 /MD /W3 /GX /O2 /I "..\..\visit\databases\PDB" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "GENERAL_PLUGIN_EXPORTS" /D "MDSERVER_PLUGIN_EXPORTS" /D "USING_MSVC6" /YX /FD /TP /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib state.lib misc.lib plugin.lib pipeline_ser.lib dbatts.lib database_ser.lib avtexceptions.lib vtkCommon.lib silohdf5.lib /nologo /dll /machine:I386 /out:"Release/libMPDBDatabase.dll"
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy Release\libMPDBDatabase.dll ..\..\bin\Release\databases
# End Special Build Tool

!ELSEIF  "$(CFG)" == "PDBM - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "SILOM_EXPORTS" /YX /FD /GZ /c
# ADD CPP /nologo /G6 /MDd /W3 /Gm /GX /ZI /Od /I "..\..\visit\databases\PDB" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "GENERAL_PLUGIN_EXPORTS" /D "MDSERVER_PLUGIN_EXPORTS" /D "USING_MSVC6" /YX /FD /GZ /TP /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib state.lib misc.lib plugin.lib pipeline_ser.lib dbatts.lib database_ser.lib avtexceptions.lib vtkCommon.lib silohdf5.lib /nologo /dll /debug /machine:I386 /out:"Debug/libMPDBDatabase.dll" /pdbtype:sept
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy Debug\libMPDBDatabase.dll ..\..\bin\Debug\databases
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "PDBM - Win32 Release"
# Name "PDBM - Win32 Debug"
# Begin Source File

SOURCE=..\..\visit\databases\PDB\avtPDBFileFormat.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\LEOSFileFormat.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\LEOSFileReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\MaterialEncoder.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\PDBCommonPluginInfo.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\PDBFileObject.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\PDBMDServerPluginInfo.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\PDBPluginInfo.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\PDBReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\PF3DFileFormat.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\PP_Z_MTSD_FileFormat.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\PP_Z_STSD_FileFormat.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\databases\PDB\PP_ZFileReader.C
# End Source File
# End Target
# End Project
