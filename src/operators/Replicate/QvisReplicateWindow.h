/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
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

#ifndef QVISREPLICATEWINDOW_H
#define QVISREPLICATEWINDOW_H

#include <QvisOperatorWindow.h>
#include <AttributeSubject.h>

class ReplicateAttributes;
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

// ****************************************************************************
// Class: QvisReplicateWindow
//
// Purpose: 
//   Defines QvisReplicateWindow class.
//
// Notes:      This class was automatically generated!

// Programmer: xml2window
// Creation:   Thu Mar 22 12:57:41 PDT 2007
//
// Modifications:
//   
// ****************************************************************************

class QvisReplicateWindow : public QvisOperatorWindow
{
    Q_OBJECT
  public:
    QvisReplicateWindow(const int type,
                         ReplicateAttributes *subj,
                         const char *caption = 0,
                         const char *shortName = 0,
                         QvisNotepadArea *notepad = 0);
    virtual ~QvisReplicateWindow();
    virtual void CreateWindowContents();
  protected:
    void UpdateWindow(bool doAll);
    virtual void GetCurrentValues(int which_widget);
  private slots:
    void useUnitCellVectorsChanged(bool val);
    void xVectorProcessText();
    void yVectorProcessText();
    void zVectorProcessText();
    void xReplicationsProcessText();
    void yReplicationsProcessText();
    void zReplicationsProcessText();
    void mergeResultsChanged(bool val);
    void replicateUnitCellAtomsChanged(bool val);
  private:
    QCheckBox *useUnitCellVectors;
    QLineEdit *xVector;
    QLineEdit *yVector;
    QLineEdit *zVector;
    QLineEdit *xReplications;
    QLineEdit *yReplications;
    QLineEdit *zReplications;
    QCheckBox *mergeResults;
    QCheckBox *replicateUnitCellAtoms;
    QLabel *useUnitCellVectorsLabel;
    QLabel *xVectorLabel;
    QLabel *yVectorLabel;
    QLabel *zVectorLabel;
    QLabel *xReplicationsLabel;
    QLabel *yReplicationsLabel;
    QLabel *zReplicationsLabel;
    QLabel *mergeResultsLabel;
    QLabel *replicateUnitCellAtomsLabel;

    ReplicateAttributes *atts;
};



#endif
