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

#ifndef QVIS_EXPRESSIONS_WINDOW_H
#define QVIS_EXPRESSIONS_WINDOW_H
#include <gui_exports.h>
#include <QvisPostableWindowObserver.h>
#include <AttributeSubject.h>

// Forward declarations
class ExpressionList;
class QVBoxLayout;
class QHBoxLayout;
class QPushButton;
class QCheckBox;
class QGroupBox;
class QLineEdit;
class QMultiLineEdit;
class QLabel;
class QListBox;
class QComboBox;
class QvisVariableButton;

// ****************************************************************************
// Class: QvisExpressionsWindow
//
// Purpose:
//   This class contains the widgets that manipulate expressions
//
// Notes:      
//
// Programmer: Jeremy Meredith
// Creation:   October 10, 2004
//
// Modifications:
//    Jeremy Meredith, Mon Oct 25 11:32:14 PDT 2004
//    Reversed the sense of the "hidden" button.  Added a list-box-index to
//    expresion-list-index map so we didn't have to index expressions by name.
//
//    Brad Whitlock, Thu Dec 9 10:15:25 PDT 2004
//    Added a newExpression slot function and a variable button that lets
//    us pick variables from the active source.
//
// ****************************************************************************

class GUI_API QvisExpressionsWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
  public:

    QvisExpressionsWindow(ExpressionList * exprAtts_,
                          const char *caption = 0,
                          const char *shortName = 0,
                          QvisNotepadArea * notepad = 0);
    virtual ~ QvisExpressionsWindow();
    virtual void CreateWindowContents();
  public  slots:
    virtual void apply();
    void    newExpression();
  protected:
    void    UpdateWindow(bool doAll);
    void    Apply(bool forceUpdate = false);
    void    BlockAllSignals(bool);
  private slots:
    void    addExpression();
    void    delExpression();
    void    nameTextChanged(const QString&);
    void    definitionTextChanged();
    void    typeChanged(int);
    void    notHiddenChanged();
    void    displayAllVarsChanged();
    void    insertFunction(int);
    void    insertVariable(const QString &);

    void    UpdateWindowSingleItem();
    void    UpdateWindowSensitivity();

  private:
    // Widgets and layouts.
    QListBox           *exprListBox;

    QLabel             *nameEditLabel;
    QLabel             *definitionEditLabel;
    QLabel             *typeLabel;

    QLineEdit          *nameEdit;
    QComboBox          *typeList;
    QCheckBox          *notHidden;
    QMultiLineEdit     *definitionEdit;

    QPushButton        *newButton;
    QPushButton        *delButton;

    QPushButton        *insertFunctionButton;
    QPopupMenu         *insertFunctionMenu;
    QvisVariableButton *insertVariableButton;

    QCheckBox          *displayAllVars;

    // State information
    ExpressionList *exprList;
    std::map<int,int> indexMap;
};

#endif
