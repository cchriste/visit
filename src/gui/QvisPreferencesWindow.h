/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef QVIS_PREFERENCES_WINDOW_H
#define QVIS_PREFERENCES_WINDOW_H
#include <gui_exports.h>
#include <QvisPostableWindowObserver.h>
#include <TimeFormat.h>

class GlobalAttributes;
class QButtonGroup;
class QCheckBox;
class QSpinBox;

// ****************************************************************************
// Class: QvisPreferencesWindow
//
// Purpose: 
//   Defines QvisPreferencesWindow class.
//
// Programmer: Eric Brugger
// Creation:   Thu Mar 13 11:13:18 PST 2003
//
// Modifications:
//   Brad Whitlock, Fri Sep 5 15:41:31 PST 2003
//   Added a toggle for post windows when shown.
//
//   Brad Whitlock, Mon Oct 13 16:52:39 PST 2003
//   Added a toggle for timeState display mode.
//
//   Brad Whitlock, Fri Jan 30 14:23:26 PST 2004
//   I added a toggle for showing the selected files list.
//
//   Brad Whitlock, Fri Apr 9 14:12:16 PST 2004
//   I added a toggle for highlighting the selected files.
//
//   Brad Whitlock, Fri Aug 6 09:20:21 PDT 2004
//   I added toggles that let you set the prompting behavior for "make default"
//   and "automatically apply operator".
//
//   Mark C. Miller, Wed Jun  1 11:12:25 PDT 2005
//   Added setTryHarderCyclesTimes check box
//
//   Mark C. Miller, Mon Jun 11 17:45:24 PDT 2007
//   Added treatAllDBsAsTimeVarying check box
// 
//   Kathleen Bonnell, Tue Oct  9 14:40:10 PDT 2007
//   Added 'Update' method, as this class now observes GlobalAttributes.
//   Added createMeshQuality, createTimeDerivative buttons.
//
//   Cyrus Harrison, Wed Nov 28 13:28:47 PST 2007
//   Added createVectorMagnitudeToggle check box
//
//   Cyrus Harrison, Wed Nov 28 13:28:47 PST 2007
//   Removed apply functions b/c the apply button was removed.
//
//   Brad Whitlock, Thu Jan 24 11:21:56 PDT 2008
//   Added newPlotsInheritSILRestriction check box.
//
//   Brad Whitlock, Thu Jan 31 10:20:52 PST 2008
//   Added session related options.
//
//   Brad Whitlock, Wed Apr  9 11:52:15 PDT 2008
//   QString for caption, shortName.
//
//   Mark C. Miller, Tue Jun 10 22:36:25 PDT 2008
//   Added preference to ignore extents.
//
//   Hank Childs, Wed Mar 17 20:13:21 PDT 2010
//   Added preference to expand new plots.
//
//   Brad Whitlock, Fri May  7 14:29:53 PDT 2010
//   I transplanted some replace plots coding.
//
//   Eric Brugger, Tue Aug 24 12:18:44 PDT 2010
//   I added a preference to enable warning message popups.
//
//   Kathleen Biagas, Wed Aug  7 13:07:12 PDT 2013
//   Added a preference for floating point precision.
//
//   David Camp, Thu Aug  8 08:50:06 PDT 2013
//   Added the restore from last session feature. 
//
//   Cameron Christensen, Tuesday, June 10, 2014
//   Added a preference for setting the backend type.
//
// ****************************************************************************

class GUI_API QvisPreferencesWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
public:
    QvisPreferencesWindow(GlobalAttributes *subj,
                          const QString &caption = QString::null,
                          const QString &shortName = QString::null,
                          QvisNotepadArea *notepad = 0);
    virtual ~QvisPreferencesWindow();
    virtual void CreateWindowContents();
    virtual void Update(Subject *);

    void SetTimeStateFormat(const TimeFormat &fmt);
    void SetShowSelectedFiles(bool val);
    void SetAllowFileSelectionChange(bool val);
    void SetEnableWarningPopups(bool val);

    bool GetEnableWarningPopups();
signals:
    void changeTimeFormat(const TimeFormat &);
    void showSelectedFiles(bool);
    void allowFileSelectionChange(bool);
    void enableWarningPopups(bool);
protected:
    void UpdateWindow(bool doAll);
private slots:
    void cloneWindowOnFirstRefToggled(bool val);
    void postWindowsWhenShownToggled(bool val);
    void makeDefaultConfirmToggled(bool val);
    void tryHarderCyclesTimesToggled(bool val);
    void automaticallyApplyOperatorToggled(bool val);
    void handleTimeStateDisplayModeChange(int val);
    void timeStateNDigitsChanged(int val);
    void selectedFilesToggled(bool);
    void allowFileSelectionChangeToggled(bool);
    void treatAllDBsAsTimeVaryingToggled(bool);
    void createMeshQualityToggled(bool);
    void createTimeDerivativeToggled(bool);
    void createVectorMagnitudeToggled(bool);
    void newPlotsInheritSILRestrictionToggled(bool);
    void expandNewPlotsToggled(bool);
    void userDirForSessionFilesToggled(bool);
    void saveCrashRecoveryFileToggled(bool);
    void ignoreDbExtentsToggled(bool val);
    void replacePlotsToggled(bool);
    void enableWarningPopupsToggled(bool);
    void userRestoreSessionFileToggled(bool);
    void precisionTypeChanged(int);
    void backendTypeChanged(int);
private:
    QCheckBox        *cloneWindowOnFirstRefToggle;
    QCheckBox        *postWindowsWhenShownToggle;
    QCheckBox        *makeDefaultConfirmToggle;
    QCheckBox        *tryHarderCyclesTimesToggle;
    QCheckBox        *automaticallyApplyOperatorToggle;
    QCheckBox        *selectedFilesToggle;
    QCheckBox        *allowFileSelectionChangeToggle;
    QCheckBox        *treatAllDBsAsTimeVaryingToggle;
    QButtonGroup     *timeStateDisplayMode;
    QSpinBox         *timeStateNDigits;
    QCheckBox        *createMeshQualityToggle;
    QCheckBox        *createTimeDerivativeToggle;
    QCheckBox        *createVectorMagnitudeToggle;
    QCheckBox        *newPlotsInheritSILRestrictionToggle;
    QCheckBox        *expandNewPlotsToggle;
    QCheckBox        *userDirForSessionFilesToggle;
    QCheckBox        *saveCrashRecoveryFileToggle;
    QCheckBox        *ignoreDbExtentsToggle;
    QCheckBox        *replacePlotsToggle;
    QCheckBox        *enableWarningPopupsToggle;
    QCheckBox        *userRestoreSessionFileToggle;
    QButtonGroup     *precisionType;
    QButtonGroup     *backendType;
    GlobalAttributes *atts;

    TimeFormat        tsFormat;
    bool              showSelFiles;
    bool              allowFileSelChange;
    bool              enableWarnPopups;
};

#endif
