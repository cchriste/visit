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

#ifndef QVIS_COMMAND_WINDOW_H
#define QVIS_COMMAND_WINDOW_H
#include <QvisPostableWindow.h>

class QButtonGroup;
class QCheckBox;
class QComboBox;
class QPushButton;
class QTabWidget;
class QTextEdit;
class QSyntaxHighlighter;
class QVBox;

// ****************************************************************************
// Class: QvisCommandWindow
//
// Purpose:
//   This class implements a window that lets you type commands to be
//   interpreted.
//
// Notes:
//
// Programmer: Brad Whitlock
// Creation:   Mon May 9 10:20:44 PDT 2005
//
// Modifications:
//   Brad Whitlock, Fri Jan 6 13:35:40 PST 2006
//   Added new buttons for recording macros.
//
//   Brad Whitlock, Fri Jun 15 13:31:37 PST 2007
//   Added a new tab for visitrc macros.
//
//   Brad Whitlock, Wed Apr  9 11:35:25 PDT 2008
//   QString for captionString, shortName.
//
//   Cyrus Harrison, Tue Jun 10 15:00:05 PDT 2008
//   Initial Qt4 Port.
//
//   Cyrus Harrison, Mon Feb  8 15:01:39 PST 2010
//   Added syntax highlighters.
//
// ****************************************************************************

class QvisCommandWindow : public QvisPostableWindow
{
    Q_OBJECT
public:
    QvisCommandWindow(const QString &captionString = QString::null,
                      const QString &shortName = QString::null,
                      QvisNotepadArea *n = 0);
    virtual ~QvisCommandWindow();
    virtual void CreateWindowContents();
    virtual void CreateNode(DataNode *);
    virtual void SetFromNode(DataNode *, const int *borders);
signals:
    void runCommand(const QString &);
public slots:
    void acceptRecordedMacro(const QString &);
private slots:
    void executeClicked(int);
    void clearClicked(int);

    void macroRecordClicked();
    void macroPauseClicked();
    void macroEndClicked();
    void macroAppendClicked(bool);
    void macroStorageActivated(int);
    void macroClearClicked();
    void macroUpdateClicked();
    void macroCreate(int);

    void textChanged0();
    void textChanged1();
    void textChanged2();
    void textChanged3();
    void textChanged4();
    void textChanged5();
    void textChanged6();
    void textChanged7();
private:
    QString fileName(int index) const;
    QString RCFileName() const;
    void LoadScripts();
    void SaveScripts();
    void UpdateMacroCheckBoxes();
    void CreateMacroFromText(const QString &);

    QTabWidget      *tabWidget;
    QButtonGroup    *executeButtonsGroup;
    QPushButton    **executeButtons;
    QButtonGroup    *clearButtonsGroup;
    QPushButton    **clearButtons;
    QButtonGroup    *addMacroButtonsGroup;
    QPushButton    **addMacroButtons;

    QTextEdit          **editors;
    QSyntaxHighlighter **highlighters;

    QPushButton     *macroRecord;
    QPushButton     *macroPause;
    QPushButton     *macroEnd;
    QCheckBox       *macroAppendCheckBox;
    QComboBox       *macroStorageComboBox;

    QWidget            *macroTab;
    QTextEdit          *macroEdit;
    QPushButton        *macroUpdateButton;
    QPushButton        *macroClearButton;
    QSyntaxHighlighter *macroHighlighter;

    int              macroStorageMode;
    bool             macroAppend;
    int              maxUserMacro;
};

#endif
