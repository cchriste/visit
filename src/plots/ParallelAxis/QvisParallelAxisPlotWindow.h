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

#ifndef QVIS_PARALLEL_AXIS_PLOT_WINDOW_H
#define QVIS_PARALLEL_AXIS_PLOT_WINDOW_H

#include <QvisPostableWindowObserver.h>
#include <ObserverToCallback.h>

#include <AttributeSubject.h>

#include <vectortypes.h>

#include <vector>
#include <string>

class ParallelAxisAttributes;
class QButtonGroup;
class QPushButton;
class QvisVariableButton;
class QLabel;


// ****************************************************************************
// Class: QvisParallelAxisPlotWindow
//
// Purpose: GUI window for the ParallelAxis (parallel coordinate) plot.
//
// Programmer: Mark Blair
// Creation:   Mon Mar 27 18:24:00 PST 2006
//
// Modifications:
//
//      Mark Blair, Wed Aug 16 17:12:00 PDT 2006
//      Removed widgets that display axis extents and extents selected by
//      Extents tool.  These were considered unnecessary.
//
// ****************************************************************************

class QvisParallelAxisPlotWindow : public QvisPostableWindowObserver
{
    Q_OBJECT
public:
    QvisParallelAxisPlotWindow(const int type,
                                ParallelAxisAttributes *parAxisAtts_,
                                const char *caption = 0,
                                const char *shortName = 0,
                                QvisNotepadArea *notepad = 0);
    virtual ~QvisParallelAxisPlotWindow();
    virtual void CreateWindowContents();

public slots:
    virtual void apply();
    virtual void makeDefault();
    virtual void reset();

protected:
    void UpdateWindow(bool doAll);
    void GetCurrentValues(int which_widget);
    void Apply(bool ignore = false);

private slots:
    void prevAxisClicked();
    void nextAxisClicked();
    void axisAdded(const QString &axisToAdd);
    void axisDeleted(const QString &axisToDelete);
    void leftAxisSelected(const QString &axisToSelect);

private:
    void UpdateShownFields(bool applyvalues);

    int                         plotType;

    ParallelAxisAttributes     *parAxisAtts;

    QLabel                      *axisVariable;
    QLabel                      *axisPosition;
    QPushButton                 *showPrevAxis;
    QPushButton                 *showNextAxis;
    QvisVariableButton          *addAxis;
    QvisVariableButton          *deleteAxis;
    QvisVariableButton          *leftAxis;

    int                          latestGUIShownOrd;
};


#endif
