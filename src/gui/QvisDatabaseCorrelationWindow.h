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

#ifndef QVIS_DATABASECORRELATION_WINDOW_H
#define QVIS_DATABASECORRELATION_WINDOW_H
#include <QvisWindowBase.h>

class DatabaseCorrelation;
class QComboBox;
class QLineEdit;
class QListWidget;
class QPushButton;

// ****************************************************************************
// Class: QvisDatabaseCorrelationWindow
//
// Purpose:
//   This class provides controls for designing a database correlation.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Sat Jan 31 02:45:19 PDT 2004
//
// Modifications:
//   Brad Whitlock, Wed Apr  9 11:48:45 PDT 2008
//   QString for caption.
//
// ****************************************************************************

class GUI_API QvisDatabaseCorrelationWindow : public QvisWindowBase
{
    Q_OBJECT
public:
    QvisDatabaseCorrelationWindow(const QString &correlationName,
        const QString &caption);
    QvisDatabaseCorrelationWindow(const DatabaseCorrelation &correlation,
        const QString &caption);
    virtual ~QvisDatabaseCorrelationWindow();

signals:
    void deleteMe(QvisWindowBase *);
protected slots:
    void setAddButtonEnabled(int);
    void setRemoveButtonEnabled(int);
    void addSources();
    void removeSources();
    void actionClicked();
    void cancelClicked();
protected:
    void CreateWidgets(const DatabaseCorrelation &correlation);
    void UpdateAddRemoveButtonsEnabledState();
    int  SelectedCount(const QListWidget *) const;
    void TransferItems(QListWidget *srcLB, QListWidget *destLB);

    static int   instanceCount;
    bool         createMode;

    // Widgets and layout
    QLineEdit    *correlationNameLineEdit;
    QListWidget     *sourcesListBox;
    QListWidget     *correlatedSourcesListBox;
    QPushButton  *addButton;
    QPushButton  *removeButton;
    QComboBox    *correlationMethodComboBox;
//    QCheckBox    *automaticNStatesCheckBox;
//    QSpinBox     *nStatesSpinBox;
};

#endif
