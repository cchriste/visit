# Microsoft Developer Studio Project File - Name="vtk_sl_io" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102

CFG=vtk_sl_io - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "vtk_sl_io.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "vtk_sl_io.mak" CFG="vtk_sl_io - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "vtk_sl_io - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "vtk_sl_io - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "vtk_sl_io - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release\vtk_sl_io"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "vtk_sl_io_EXPORTS" /YX /FD /c
# ADD CPP /nologo /G6 /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "USING_MSVC6" /D "VTK_SL_IO_EXPORTS" /YX /FD /TP /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib vtkCommon.lib vtkFiltering.lib vtkIO.lib vtkGraphics.lib vtkHybrid.lib vtkRendering.lib misc.lib avtexceptions.lib vtktiff.lib /nologo /dll /machine:I386
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy Release\vtk_sl_io.lib ..\..\lib\Release\vtk_sl_io.lib	copy Release\vtk_sl_io.dll ..\..\bin\Release\vtk_sl_io.dll
# End Special Build Tool

!ELSEIF  "$(CFG)" == "vtk_sl_io - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug\vtk_sl_io"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "vtk_sl_io_EXPORTS" /YX /FD /GZ /c
# ADD CPP /nologo /G6 /MDd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "USING_MSVC6" /D "VTK_SL_IO_EXPORTS" /YX /FD /GZ /TP /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib vtkCommon.lib vtkFiltering.lib vtkIO.lib vtkGraphics.lib vtkHybrid.lib vtkRendering.lib misc.lib avtexceptions.lib vtktiff.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Cmds=copy Debug\vtk_sl_io.dll ..\..\bin\Debug\vtk_sl_io.dll	copy Debug\vtk_sl_io.lib ..\..\lib\Debug\vtk_sl_io.lib
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "vtk_sl_io - Win32 Release"
# Name "vtk_sl_io - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkAppendFilter.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkAppendPolyData.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkDataObjectReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkDataObjectWriter.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkDataReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkDataSetReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkDataSetWriter.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkDataWriter.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkOBJReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkPolyDataReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkPolyDataWriter.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkRectilinearGridReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkRectilinearGridWriter.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkSTLReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkStructuredGridReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkStructuredGridWriter.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkStructuredPointsReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkStructuredPointsWriter.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkUnstructuredGridReader.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkUnstructuredGridWriter.C
# End Source File
# Begin Source File

SOURCE=..\..\visit\visit_vtk\sl_io\vtkWriter.C
# End Source File
# End Group
# End Target
# End Project
