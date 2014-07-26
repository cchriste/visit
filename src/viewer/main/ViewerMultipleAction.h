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

#ifndef VIEWER_MULTIPLE_ACTION_H
#define VIEWER_MULTIPLE_ACTION_H
#include <viewer_exports.h>
#include <ViewerActionBase.h>

#include <QIcon>
#include <QPixmap>
#include <QString>

#include <vector>

class QAction;
class QActionGroup;
class QMenu;
class QToolBar;

// ****************************************************************************
// Class: ViewerMultipleAction
//
// Purpose:
//   This is a base class for actions that multiple toolbar buttons or menu
//   options but still service a single RPC.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Jan 30 13:49:49 PST 2003
//
// Modifications:
//   Sean Ahern, Thu Feb 20 01:51:35 America/Los_Angeles 2003
//   Added the ability to set large and small icons.
//
//   Brad Whitlock, Tue Aug 26 17:12:02 PST 2003
//   I added the isExclusive flag.
//
//   Brad Whitlock, Tue Apr 29 11:16:17 PDT 2008
//   Converted to QString for menu items.
//
//   Brad Whitlock, Thu May 22 13:44:11 PDT 2008
//   Qt 4.
//
// ****************************************************************************

class VIEWER_API ViewerMultipleAction : public ViewerActionBase
{
    Q_OBJECT

    typedef std::vector<QAction *> ActionPointerVector;
public:
    ViewerMultipleAction(ViewerWindow *win);
    virtual ~ViewerMultipleAction();
    
    virtual void Setup();
    virtual void Execute();
    virtual void Update();

    virtual void Execute(int) = 0;

    virtual bool Enabled() const;
    virtual bool ChoiceEnabled(int i) const;
    virtual bool ChoiceChecked(int i) const;

    virtual void HideChoice(int i); 

    // Methods to add the action to the menu and toolbar.
    virtual void ConstructMenu(QMenu *menu);
    virtual void RemoveFromMenu(QMenu *menu);
    virtual void ConstructToolbar(QToolBar *toolbar);
    virtual void RemoveFromToolbar(QToolBar *toolbar);
    virtual void UpdateConstruction();

    // Methods to set the action's attributes.
    virtual void SetAllText(const QString &text);
    virtual void SetText(const QString &text);
    virtual void SetMenuText(const QString &text);
    virtual void SetToolTip(const QString &text);
    virtual void SetIcon(const QIcon &icon);

    virtual void AddChoice(const QString &menuText);
    virtual void AddChoice(const QString &menuText, const QString &toolTip, const
                           QPixmap &icon);
    virtual void AddChoice(const QString &menuText, const QString &toolTip,
                           const QPixmap &small_icon,
                           const QPixmap &large_icon);

    virtual void SetExclusive(bool val);
protected slots:
    void ActivateHelper(QAction*);
protected:
    bool                 iconSpecified;
    int                  activeAction;
    bool                 toggled;
    QString              text;
    QString              menuText;
    QString              toolTip;
    QIcon                icon;
    QActionGroup        *action;
    QMenu               *actionMenu;
    bool                 isExclusive;
    ActionPointerVector  children;
};

#endif
