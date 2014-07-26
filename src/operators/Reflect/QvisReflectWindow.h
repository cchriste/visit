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

#ifndef QVISREFLECTWINDOW_H
#define QVISREFLECTWINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>

class ReflectAttributes;
class QLabel;
class QCheckBox;
class QLineEdit;
class QVBox;
class QButtonGroup;
class QComboBox;
class QRadioButton;
class QvisReflectWidget;

// ****************************************************************************
//  Class:  QvisReflectWindow
//
//  Purpose:
//    The attribute window for the Reflect operator
//
//  Programmer:  Jeremy Meredith
//  Creation:    March 11, 2002
//
//  Modifications:
//    Brad Whitlock, Fri Apr 12 13:12:17 PST 2002
//    Made it inherit from QvisOperatorWindow.
//
//    Brad Whitlock, Mon Mar 3 11:40:37 PDT 2003
//    I spruced up the window.
//
//    Brad Whitlock, Wed Jun 25 09:22:58 PDT 2003
//    I added a 2D view of the window.
//
// ****************************************************************************

class QvisReflectWindow : public QvisOperatorWindow
{
    Q_OBJECT
  public:
    QvisReflectWindow(const int type,
                         ReflectAttributes *subj,
                         const QString &caption = QString::null,
                         const QString &shortName = QString::null,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisReflectWindow();
    virtual void CreateWindowContents();
  protected:
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
    void UpdateOctantMenuContents();
  private slots:
    void octantChanged(int val);
    void xBoundaryChanged(int val);
    void yBoundaryChanged(int val);
    void zBoundaryChanged(int val);
    void specifiedXProcessText();
    void specifiedYProcessText();
    void specifiedZProcessText();
    void selectOctants(bool *octants);
    void selectMode(int);
  private:
    bool               userSetMode;
    bool               mode2D;

    QButtonGroup      *modeButtons;
    QLabel            *originalDataLabel;
    QComboBox         *octant;
    QLabel            *reflectionLabel;
    QvisReflectWidget *reflect;
    QButtonGroup      *xBound;
    QRadioButton      *xUseData;
    QRadioButton      *xSpecify;
    QLineEdit         *specifiedX;
    QButtonGroup      *yBound;
    QRadioButton      *yUseData;
    QRadioButton      *ySpecify;
    QLineEdit         *specifiedY;
    QButtonGroup      *zBound;
    QRadioButton      *zUseData;
    QRadioButton      *zSpecify;
    QLineEdit         *specifiedZ;

    ReflectAttributes *atts;
};



#endif
