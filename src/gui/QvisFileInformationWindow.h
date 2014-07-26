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

#ifndef QVIS_FILE_INFORMATION_WINDOW_H
#define QVIS_FILE_INFORMATION_WINDOW_H
#include <gui_exports.h>
#include <QvisPostableWindowObserver.h>

class FileServerList;
class QTextEdit;

// ****************************************************************************
// Class: QvisFileInformationWindow
//
// Purpose:
//   Observes the file server and updates its contents when the open file
//   changes.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Aug 20 13:51:40 PST 2002
//
// Modifications:
//   Brad Whitlock, Wed Apr  9 10:46:33 PDT 2008
//   QString for caption and shortName.
//
// ****************************************************************************

class GUI_API QvisFileInformationWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
public:
    QvisFileInformationWindow(FileServerList *fs, const QString &caption = QString::null,
        const QString &shortName = QString::null, QvisNotepadArea *notepad = 0);
    virtual ~QvisFileInformationWindow();
    virtual void CreateWindowContents();
protected:
    virtual void UpdateWindow(bool doAll);
private:
    static const char *titleSeparator;
    QTextEdit    *outputText;
};

#endif
