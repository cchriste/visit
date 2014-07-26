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

#ifndef SPLASHSCREEN_H
#define SPLASHSCREEN_H

#include <QFrame>
#include <vector>

// Forward declarations.
class QLabel;
class QPushButton;
class QProgressBar;
class QPainter;
class QTimer;
class QVBoxLayout;

// ****************************************************************************
// Class: SplashScreen
//
// Purpose:
//   This is class creates a splashscreen.
//
// Notes:      
//
// Programmer: Sean Ahern
// Creation:   Thu Sep 6 16:34:45 PST 2001
//
// Modifications:
//   Brad Whitlock, Thu Sep 6 16:35:18 PST 2001
//   Modified the widget so it observes two subjects.
//
//   Brad Whitlock, Wed Jun 18 17:53:09 PST 2003
//   I changed this class so it is a simple widget.
//
//   Brad Whitlock, Tue Jan  8 13:43:02 PST 2008
//   Added signals for showing copyright and contributors.
//
//   Brad Whitlock, Wed Apr  9 10:26:06 PDT 2008
//   Use QString instead of const char *.
//
//   Brad Whitlock, Fri May 30 15:21:37 PDT 2008
//   Qt 4.
//
// ****************************************************************************

class SplashScreen : public QFrame
{
    Q_OBJECT
public:
    SplashScreen(bool cyclePictures = false);
    ~SplashScreen();

    void Progress(const QString &msg, int progress);
    void SetDisplayAsSplashScreen(bool mode);
    void About();
signals:
    void showCopyright();
    void showContributors();
public slots:
    virtual void show();
    virtual void hide();
    void nextPicture();
private slots:
    void emitShowCopyright();
    void emitShowContributors();
protected:
    void CreateAboutButtons();

    QLabel                 *pictureLabel;
    std::vector<QPixmap>   pictures;
    int                    curPicture;
    bool                   splashMode;
    QTimer                 *timer;
    QLabel                 *text;
    QProgressBar           *progress;
    QPushButton            *dismissButton;
    QPushButton            *copyrightButton;
    QPushButton            *contributorButton;
    QVBoxLayout            *topLayout;
    QVBoxLayout            *lLayout;
    QVBoxLayout            *rLayout;
};

#endif
