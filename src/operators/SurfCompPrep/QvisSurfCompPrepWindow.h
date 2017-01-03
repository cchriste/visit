/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
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

#ifndef QVISSURFCOMPPREPWINDOW_H
#define QVISSURFCOMPPREPWINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>

class SurfCompPrepAttributes;
class QLabel;
class QLineEdit;
class QButtonGroup;

// ****************************************************************************
// Class: QvisSurfCompPrepWindow
//
// Purpose:
//    Defines QvisSurfCompPrepWindow class.
//
// Notes:      Autogenerated by xml2window.
//
// Programmer: xml2window
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class QvisSurfCompPrepWindow : public QvisOperatorWindow
{
    Q_OBJECT
  public:
    QvisSurfCompPrepWindow(const int type,
                         SurfCompPrepAttributes *subj,
                         const QString &caption = QString::null,
                         const QString &shortName = QString::null,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisSurfCompPrepWindow();
    virtual void CreateWindowContents();
  protected:
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
  private slots:
    void surfaceTypeChanged(int val);
    void coordSystemChanged(int val);
    void thetaStartProcessText();
    void thetaStopProcessText();
    void thetaStepsProcessText();
    void phiStartProcessText();
    void phiStopProcessText();
    void phiStepsProcessText();
    void startRadiusProcessText();
    void endRadiusProcessText();
    void radiusStepsProcessText();
    void xStartProcessText();
    void xStopProcessText();
    void xStepsProcessText();
    void yStartProcessText();
    void yStopProcessText();
    void yStepsProcessText();
    void zStartProcessText();
    void zStopProcessText();
    void zStepsProcessText();
  private:
    QButtonGroup *surfaceType;
    QButtonGroup *coordSystem;
    QLineEdit *thetaStart;
    QLineEdit *thetaStop;
    QLineEdit *thetaSteps;
    QLineEdit *phiStart;
    QLineEdit *phiStop;
    QLineEdit *phiSteps;
    QLineEdit *startRadius;
    QLineEdit *endRadius;
    QLineEdit *radiusSteps;
    QLineEdit *xStart;
    QLineEdit *xStop;
    QLineEdit *xSteps;
    QLineEdit *yStart;
    QLineEdit *yStop;
    QLineEdit *ySteps;
    QLineEdit *zStart;
    QLineEdit *zStop;
    QLineEdit *zSteps;
    QLabel *surfaceTypeLabel;
    QLabel *coordSystemLabel;
    QLabel *thetaStartLabel;
    QLabel *thetaStopLabel;
    QLabel *thetaStepsLabel;
    QLabel *phiStartLabel;
    QLabel *phiStopLabel;
    QLabel *phiStepsLabel;
    QLabel *startRadiusLabel;
    QLabel *endRadiusLabel;
    QLabel *radiusStepsLabel;
    QLabel *xStartLabel;
    QLabel *xStopLabel;
    QLabel *xStepsLabel;
    QLabel *yStartLabel;
    QLabel *yStopLabel;
    QLabel *yStepsLabel;
    QLabel *zStartLabel;
    QLabel *zStopLabel;
    QLabel *zStepsLabel;

    SurfCompPrepAttributes *atts;
};



#endif
