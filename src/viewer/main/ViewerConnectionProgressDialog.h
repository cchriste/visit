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

#ifndef VIEWER_CONNECTION_PROGRESS_DIALOG_H
#define VIEWER_CONNECTION_PROGRESS_DIALOG_H
#include <QWidget>

// Forward declarations.
class QLabel;
class QPushButton;
class QRadioButton;
class QTimer;

// ****************************************************************************
// Class: ViewerViewerConnectionProgressDialog
//
// Purpose:
//   This class contains a progress dialog that is displayed when we launch
//   processes and it takes a long time.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Sep 26 13:20:25 PST 2002
//
// Modifications:
//   Brad Whitlock, Thu Oct 24 11:46:47 PDT 2002
//   I increased the default timeout to 6 seconds.
//
//   Brad Whitlock, Wed May 7 10:14:04 PDT 2003
//   I added methods to set/get the component name. I also added methods
//   to set/get a flag that determines if the window should be hidden.
//
//   Brad Whitlock, Tue Nov 29 16:01:20 PST 2011
//   I added various methods for getting/setting values and I simplified the
//   constructor.
//
// ****************************************************************************

class ViewerConnectionProgressDialog : public QWidget
{
    Q_OBJECT
public:
    ViewerConnectionProgressDialog(const QString &host);
    virtual ~ViewerConnectionProgressDialog();

    bool getCancelled() const { return cancelled; };
    void setTimeout(int t) { timeout = t; };

    const QString &getHostName() const { return hostName; };
    void setHostName(const QString &host);

    const QString &getComponentName() const { return componentName; };
    void setComponentName(const QString &cn);

    void setParallel(bool);
    bool getParallel() const { return parallel; }

    bool getIgnoreHide() const   { return ignoreHide; }
    void setIgnoreHide(bool val) { ignoreHide = val; }
public slots:
    virtual void show();
    virtual void hide();
private slots:
    void updateAnimation();
    void timedShow();
    void cancel();
private:
    void updateMessage();

    QTimer       *timer;
    QPushButton  *cancelButton;
    QRadioButton *rb[6];
    QLabel       *msgLabel;
    QLabel       *rightComputer;
    QString       componentName;
    QString       hostName;
    bool          parallel;
    int           iconFrame;
    bool          cancelled;
    int           timeout;
    bool          cancelledShow;
    bool          ignoreHide;
};

#endif
