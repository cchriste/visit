/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
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

#ifndef QVISCREATEBONDSWINDOW_H
#define QVISCREATEBONDSWINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>

class CreateBondsAttributes;
class QLabel;
class QCheckBox;
class QLineEdit;
class QSpinBox;
class QVBox;
class QButtonGroup;
class QvisColorTableButton;
class QvisOpacitySlider;
class QvisColorButton;
class QvisLineStyleWidget;
class QvisLineWidthWidget;
class QvisVariableButton;
class QPushButton;
class QListBox;
class QComboBox;
class QListView;
class QListViewItem;
class QvisElementButton;

// ****************************************************************************
// Class: QvisCreateBondsWindow
//
// Purpose: 
//   Defines QvisCreateBondsWindow class.
//
// Notes:      This class was automatically generated!

// Programmer: xml2window
// Creation:   Wed Apr 19 17:23:18 PST 2006
//
// Modifications:
//    Jeremy Meredith, Mon Feb 11 17:03:15 EST 2008
//    Added support for wildcards for atomic numbers.  Necessitated
//    adding an up/down button (since order is now significant).
//   
// ****************************************************************************

class QvisCreateBondsWindow : public QvisOperatorWindow
{
    Q_OBJECT
  public:
    QvisCreateBondsWindow(const int type,
                         CreateBondsAttributes *subj,
                         const char *caption = 0,
                         const char *shortName = 0,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisCreateBondsWindow();
    virtual void CreateWindowContents();
  protected:
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
    int  GetItemIndex(QListViewItem*);
    int  GetListLength();
  private slots:
    void UpdateWindowSingleItem();
    void elementVariableChanged(const QString &varName);
    void bondsListNew();
    void bondsListDel();
    void bondsListUp();
    void bondsListDown();
    void minDistTextChanged(const QString&);
    void maxDistTextChanged(const QString&);
    void minDistReturnPressed();
    void maxDistReturnPressed();
    void maxBondsReturnPressed();
    void firstElementChanged(int);
    void secondElementChanged(int);
  private:
    QvisVariableButton *elementVariable;
    QLabel *elementVariableLabel;

    QLineEdit *maxBonds;
    QLabel    *maxBondsLabel;

    QPushButton     *newButton;
    QPushButton     *delButton;

    QPushButton     *upButton;
    QPushButton     *downButton;

    QListView       *bondsList;
    CreateBondsAttributes *atts;

    QvisElementButton *firstElement;
    QvisElementButton *secondElement;
    QLineEdit         *minDist;
    QLineEdit         *maxDist;
};



#endif
