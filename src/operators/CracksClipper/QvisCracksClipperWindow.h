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

#ifndef QVISCRACKSCLIPPERWINDOW_H
#define QVISCRACKSCLIPPERWINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>

class CracksClipperAttributes;
class QCheckBox;
class QComboBox;
class QLabel;
class QvisVariableButton;

// ****************************************************************************
//  Class: QvisCracksClipperWindow
//
//  Purpose: 
//    Defines QvisCracksClipperWindow class.
//
//  Notes:      This class was automatically generated!
//
//  Programmer: xml2window
//  Creation:   Mon Aug 22 09:10:02 PDT 2005
//
//  Modifications:
//    Kathleen Bonnell, Mon May  7 15:48:42 PDT 2007
//    Added calculateDensity, inMassVar, outDenVar.
//   
//    Kathleen Bonnell, Wed Sep 29 09:01:48 PDT 2010
//    Removed calculateDensity, outDenVar.
//
// ****************************************************************************

class QvisCracksClipperWindow : public QvisOperatorWindow
{
    Q_OBJECT
  public:
    QvisCracksClipperWindow(const int type,
                         CracksClipperAttributes *subj,
                         const QString &caption = QString::null,
                         const QString &shortName = QString::null,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisCracksClipperWindow();
    virtual void CreateWindowContents();
  protected:
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
  private slots:
    void crack1VarChanged(const QString &);
    void crack2VarChanged(const QString &);
    void crack3VarChanged(const QString &);
    void strainVarChanged(const QString &);
    void showCrack1Changed(bool val);
    void showCrack2Changed(bool val);
    void showCrack3Changed(bool val);
    void inMassVarChanged(const QString &);
  private:
    QvisVariableButton *crack1Var;
    QvisVariableButton *crack2Var;
    QvisVariableButton *crack3Var;
    QvisVariableButton *strainVar;
    QCheckBox *showCrack1;
    QCheckBox *showCrack2;
    QCheckBox *showCrack3;
    QvisVariableButton *inMassVar;
    QLabel *inMassVarLabel;

    CracksClipperAttributes *atts;
};



#endif
