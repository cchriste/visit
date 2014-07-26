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

#ifndef SIMCOMMANDSLOTS_H
#define SIMCOMMANDSLOTS_H

#include <QDialog>
#include <QDateTime>
#include <QDateTimeEdit>

#include <EngineList.h>
#include <StatusAttributes.h>
#include <ViewerProxy.h>

class QEvent;

// ****************************************************************************
// Class: SimCommandSlots
//
// Purpose:
//   This class implements the slots defined for the qt ui that is loaded
//   dynamically.
//
// Notes:      
//
// Programmer: Shelly Prevost
// Creation:   Jan 6, 2006
//
// Modifications:
//
// ****************************************************************************

class SimCommandSlots:  public QObject
{
    Q_OBJECT
public:
    SimCommandSlots(ViewerProxy *theViewer, EngineList *engineList,
                    int index );
    void SetParent(QDialog * parent ) {p = parent;};

public slots:
    void SimCommands() {printf("command Slot Fired\n");};
    void ClickedHandler();
    void TextChangedHandler(const QString &);
    void ValueChangedHandler(int);
    void ValueChangedHandler(const QTime &theTime);
    void ValueChangedHandler(const QDate &theDate);
    void ActivatedHandler(int);
    void StateChangedHandler(int);
    void CurrentChangedHandler(int, int);
    void ValueChangedHandler(int, int);
    int  SendCMD(QString sig, const  QObject *ui,QString value);

private:
    bool eventFilter(QObject *, QEvent *) {printf("inside eventFilter\n");return false;};

    QDialog     *p;
    EngineList  *engines;
    ViewerProxy *viewer;
    int          simIndex;
};

#endif
