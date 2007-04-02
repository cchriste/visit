/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef QVIS_RECENT_PATH_REMOVAL_WINDOW_H
#define QVIS_RECENT_PATH_REMOVAL_WINDOW_H
#include <gui_exports.h>
#include <QvisDelayedWindowObserver.h>
#include <QualifiedFilename.h>

// Forward declarations
class QGroupBox;
class QListBox;
class QPushButton;

// ****************************************************************************
// Class: QvisRecentPathRemovalWindow
//
// Purpose:
//   Removes paths from the recently visited path list in the file selection
//   window.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Oct 10 14:57:48 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class GUI_API QvisRecentPathRemovalWindow : public QvisDelayedWindowObserver
{
    Q_OBJECT
public:
    QvisRecentPathRemovalWindow(Subject *s, const char *captionString);
    virtual ~QvisRecentPathRemovalWindow();
public slots:
    virtual void show();
protected:
    virtual void CreateWindowContents();
    virtual void UpdateWindow(bool doAll);
    void UpdateWidgets();
private slots:
    void removePaths();
    void removeAllPaths();
    void invertSelection();
    void applyDismiss();
    void handleCancel();
private:
    QGroupBox   *removalControlsGroup;
    QListBox    *removalListBox;
    QPushButton *removeButton;
    QPushButton *removeAllButton;
    QPushButton *invertSelectionButton;

    QPushButton *okButton;
    QPushButton *cancelButton;

    QualifiedFilenameVector paths;
};

#endif
