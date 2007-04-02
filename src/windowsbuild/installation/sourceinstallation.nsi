; Script generated by the HM NIS Edit Script Wizard.

##############################################################################
#
# File: sourceinstallation.nsi
#
# Purpose: This file contains the instructions that the NSIS installer needs
#          in order to create an installation program for VisIt's source code.
#
# Programmer: Brad Whitlock
# Date:       Mon, Dec 15 11:26:34 PDT 2003
#
# Modifications:
#   Brad Whitlock, Mon Feb 9 14:45:04 PST 2004
#   Updated for 1.2.7
#
#   Brad Whitlock, Fri Mar 5 09:40:50 PDT 2004
#   Updated for 1.2.8.
#
#   Brad Whitlock, Wed Apr 14 16:23:17 PST 2004
#   Updated for 1.3.
#
#   Brad Whitlock, Thu May 27 18:19:19 PST 2004
#   Updated for 1.3.1
#
#   Brad Whitlock, Tue Jun 29 13:08:33 PST 2004
#   Updated for 1.3.2.
#
#   Brad Whitlock, Wed Jul 14 09:34:28 PDT 2004
#   Updated for 1.3.3 and changed registry keys.
#
#   Brad Whitlock, Tue Aug 24 11:37:40 PDT 2004
#   Updated for 1.3.4. We now have MSVC7.Net specific stuff to distribute.
#
#   Brad Whitlock, Thu Sep 23 09:39:12 PDT 2004
#   Updated for 1.3.5
#
#   Brad Whitlock, Wed Nov 3 14:05:25 PST 2004
#   Updated for 1.4.
#
#   Brad Whitlock, Wed Jan 5 17:24:33 PST 2005
#   Updated for 1.4.1.
#
#   Brad Whitlock, Thu Feb 24 16:09:50 PST 2005
#   Updated for 1.4.2.
#
#   Brad Whitlock, Tue May 10 14:12:58 PST 2005
#   Updated for 1.4.3.
#
#   Brad Whitlock, Tue Jul 12 15:44:27 PST 2005
#   Updated for 1.4.4.
#
##############################################################################

; HM NIS Edit Wizard helper defines
!define PRODUCT_NAME "VisIt for Windows source code"
!define PRODUCT_VERSION "1.4.4"
!define PRODUCT_PUBLISHER "LLNL"
!define PRODUCT_WEB_SITE "http://www.llnl.gov/visit"
!define PRODUCT_DIR_REGKEY "Software\Microsoft\Windows\CurrentVersion\App Paths\visitdev${PRODUCT_VERSION}"

SetCompressor bzip2

; MUI 1.67 compatible ------
!include "MUI.nsh"

; MUI Settings
!define MUI_ABORTWARNING
!define MUI_ICON "..\resources\visit.ico"

; Welcome page
!insertmacro MUI_PAGE_WELCOME
; License page
!insertmacro MUI_PAGE_LICENSE "copyright.txt"
; Directory page
!insertmacro MUI_PAGE_DIRECTORY
; Instfiles page
!insertmacro MUI_PAGE_INSTFILES
; Finish page
!define MUI_FINISHPAGE_SHOWREADME "$INSTDIR\VisItBuildInstructionsOnWindows.doc"
!insertmacro MUI_PAGE_FINISH

; Language files
!insertmacro MUI_LANGUAGE "English"

; Reserve files
!insertmacro MUI_RESERVEFILE_INSTALLOPTIONS

; MUI end ------

Name "${PRODUCT_NAME} ${PRODUCT_VERSION}"
OutFile "..\installation\visitdev${PRODUCT_VERSION}.exe"
InstallDir "C:\VisItDev${PRODUCT_VERSION}"
InstallDirRegKey HKLM "${PRODUCT_DIR_REGKEY}" ""
ShowInstDetails show

Section "ProjectFiles" SEC02
  SetOutPath "$INSTDIR"
  File /r "..\projects"
  File /r "..\projects-MSVC7.Net"
SectionEnd

Section "IncludeFiles" SEC03
  SetOutPath "$INSTDIR"
  File /r "..\include"
SectionEnd

Section "InstallationFiles" SEC04
  # Get the instructions file
  SetOutPath "$INSTDIR"
  File "..\VisItBuildInstructionsOnWindows.doc"
  File "..\BUILD_NOTES.txt"
  # The installation files that we use to build the binary and source distributions.
  SetOutPath "$INSTDIR\installation"
  File "binaryinstallation.nsi"
  File "binaryinstallation-MSVC7.Net.nsi"
  File "sourceinstallation.nsi"
  File "copyright.txt"
  File "NetworkConfig.ini"
SectionEnd

Section "ScriptFiles" SEC05
  # The script files that we need to build on Windows.
  SetOutPath "$INSTDIR"
  File /r "..\script"
SectionEnd

Section "BinFiles" SEC06
  SetOutPath "$INSTDIR"
  File /r "..\bin"
SectionEnd

Section "VisItSource" SEC07
  SetOutPath "$INSTDIR"
  File /r "..\visit"
SectionEnd

Section "LibFiles" SEC08
  SetOutPath "$INSTDIR"
  File /r "..\lib"
SectionEnd

Section "Resources" SEC09
  SetOutPath "$INSTDIR"
  File /r "..\resources"
SectionEnd

Section -Post
  WriteRegStr HKLM "${PRODUCT_DIR_REGKEY}" "" "$INSTDIR\makensis.exe"
  
  # Set the VISITDEVDIR key in the registry.
  WriteRegStr HKCR "VISIT${PRODUCT_VERSION}" "VISITDEVDIR" "$INSTDIR"
  WriteRegStr HKCU "VISIT${PRODUCT_VERSION}" "VISITDEVDIR" "$INSTDIR"

  # Set the VISITDEVDIR key in the registry so it will be set as an 
  # environment variable for the current user.
  WriteRegStr HKCU "Environment" "VISITDEVDIR" "$INSTDIR"
SectionEnd
