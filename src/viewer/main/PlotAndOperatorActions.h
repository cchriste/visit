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

#ifndef PLOT_AND_OPERATOR_ACTIONS_H
#define PLOT_AND_OPERATOR_ACTIONS_H
#include <ViewerAction.h>
#include <ViewerMultipleAction.h>
#include <VariableMenuPopulator.h>
#include <vectortypes.h>

// ****************************************************************************
// Class: AddOperatorAction
//
// Purpose:
//   This action adds an operator to the plots in a window's plot list.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Mon Mar 17 09:24:01 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class AddOperatorAction : public ViewerMultipleAction
{
public:
    AddOperatorAction(ViewerWindow *win);
    virtual ~AddOperatorAction();

    virtual void Setup();

    virtual void Execute(int);
    virtual bool Enabled() const;
    virtual bool ChoiceEnabled(int i) const;

    virtual void ConstructToolbar(QToolBar *toolbar);
private:
    intVector graphicalPlugins;
};

// ****************************************************************************
// Class: PromoteOperatorAction
//
// Purpose:
//   This action promotes an operator in the window's selected plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Apr 10 09:38:34 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class PromoteOperatorAction : public ViewerAction
{
public:
    PromoteOperatorAction(ViewerWindow *);
    virtual ~PromoteOperatorAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: DemoteOperatorAction
//
// Purpose:
//   This action demotes an operator in the window's selected plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Apr 10 09:38:34 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class DemoteOperatorAction : public ViewerAction
{
public:
    DemoteOperatorAction(ViewerWindow *);
    virtual ~DemoteOperatorAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: RemoveOperatorAction
//
// Purpose:
//   This action removes an operator from the window's selected plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Apr 10 09:38:34 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class RemoveOperatorAction : public ViewerAction
{
public:
    RemoveOperatorAction(ViewerWindow *);
    virtual ~RemoveOperatorAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: RemoveLastOperatorAction
//
// Purpose:
//   This action removes the last operator from the window's selected plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Mon Mar 17 09:51:18 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class RemoveLastOperatorAction : public ViewerAction
{
public:
    RemoveLastOperatorAction(ViewerWindow *);
    virtual ~RemoveLastOperatorAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: RemoveAllOperatorsAction
//
// Purpose:
//   This action removes all operators from the window's selected plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Mon Mar 17 09:51:41 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class RemoveAllOperatorsAction : public ViewerAction
{
public:
    RemoveAllOperatorsAction(ViewerWindow *);
    virtual ~RemoveAllOperatorsAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: SetOperatorOptionsAction
//
// Purpose:
//   This action sets the operator attributes for an operator.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:48:55 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class SetOperatorOptionsAction : public ViewerAction
{
public:
    SetOperatorOptionsAction(ViewerWindow *);
    virtual ~SetOperatorOptionsAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: AddPlotAction
//
// Purpose:
//   This action adds a plot to the window's plot list.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Mon Mar 17 09:24:19 PDT 2003
//
// Modifications:
//   Brad Whitlock, Mon Sep 29 17:38:36 PST 2003
//   I separated host and database into two fields.
//
//   Brad Whitlock, Tue Mar 16 15:35:44 PST 2004
//   I added changeMenuIconSize.
//
//   Brad Whitlock, Fri Apr 15 13:57:12 PST 2005
//   I removed host and database.
//
//   Brad Whitlock, Thu May 29 16:42:43 PDT 2008
//   Removed menu and a slot.
//
//   Brad Whitlock, Fri Nov 19 15:05:34 PST 2010
//   I added CreatePlotMenu and DeletePlotMenu.
//
// ****************************************************************************

class AddPlotAction : public ViewerMultipleAction
{
    Q_OBJECT

    typedef struct
    {
        int                    index;
        QvisVariablePopupMenu *varMenu;
        int                    varTypes;
    } PluginEntry;

    typedef std::vector<PluginEntry> PluginEntryVector;
public:
    AddPlotAction(ViewerWindow *win);
    virtual ~AddPlotAction();

    virtual void Update();
    virtual void Execute(int);
    virtual bool Enabled() const;
    virtual bool ChoiceEnabled(int i) const;

    virtual void ConstructMenu(QMenu *);
    virtual void RemoveFromMenu(QMenu *);
    virtual void ConstructToolbar(QToolBar *toolbar);
private slots:
    void addPlot(int, const QString &);
    void changeMenuIconSize(bool);
private:
    void CreatePlotMenu(int);
    void DeletePlotMenu(int);
    int                   maxPixmapWidth, maxPixmapHeight;
    PluginEntryVector     pluginEntries;
    VariableMenuPopulator menuPopulator;
};

// ****************************************************************************
// Class: AddEmbeddedPlotAction
//
// Purpose:
//   This action adds a plot to the window's plot list, using an id specified
//   by the embedding code.
//   It also differs from the regular AddPlot because it isn't exposed in a menu
//
// Notes:      
//
// Programmer: Marc Durant
// Creation:   June 19, 2011
//
// Modifications:
//
// ****************************************************************************

class AddEmbeddedPlotAction : public ViewerAction
{
public:
  AddEmbeddedPlotAction(ViewerWindow *win);
  virtual ~AddEmbeddedPlotAction();  
  virtual void Execute();
  virtual bool Enabled() const;
};

// ****************************************************************************
// Class: DrawPlotsAction
//
// Purpose:
//   This action draws the window's plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Mon Mar 17 09:51:41 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class DrawPlotsAction : public ViewerAction
{
public:
    DrawPlotsAction(ViewerWindow *);
    virtual ~DrawPlotsAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: HideActivePlotsAction
//
// Purpose:
//   This action hides the window's active plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Mon Mar 17 09:51:41 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class HideActivePlotsAction : public ViewerAction
{
public:
    HideActivePlotsAction(ViewerWindow *);
    virtual ~HideActivePlotsAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: SetPlotFollowsTimeAction
//
// Purpose:
//   This action disconnects the window's active plots from the time slider.
//
// Notes:      
//
// Programmer: Ellen Tarwater
// Creation:   Thurs, Dec 6, 2007
//
// Modifications:
//   
// ****************************************************************************

class SetPlotFollowsTimeAction : public ViewerAction
{
public:
    SetPlotFollowsTimeAction(ViewerWindow *);
    virtual ~SetPlotFollowsTimeAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: DeleteActivePlotsAction
//
// Purpose:
//   This action deletes the window's active plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Mon Mar 17 09:51:41 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class DeleteActivePlotsAction : public ViewerAction
{
public:
    DeleteActivePlotsAction(ViewerWindow *);
    virtual ~DeleteActivePlotsAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: SetActivePlotsAction
//
// Purpose:
//   This action sets the window's active plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:30:16 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class SetActivePlotsAction : public ViewerAction
{
public:
    SetActivePlotsAction(ViewerWindow *);
    virtual ~SetActivePlotsAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: ChangeActivePlotsVarAction
//
// Purpose:
//   This action sets the plotted variable for the selected plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:30:16 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class ChangeActivePlotsVarAction : public ViewerAction
{
public:
    ChangeActivePlotsVarAction(ViewerWindow *);
    virtual ~ChangeActivePlotsVarAction();
    virtual void Execute();
};


// ****************************************************************************
// Class: SetPlotSILRestrictionAction
//
// Purpose:
//   This action sets the SIL restriction for a plot.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:47:41 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class SetPlotSILRestrictionAction : public ViewerAction
{
public:
    SetPlotSILRestrictionAction(ViewerWindow *);
    virtual ~SetPlotSILRestrictionAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: SetPlotOptionsAction
//
// Purpose:
//   This action sets the plot attributes for the selected plots.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:48:57 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class SetPlotOptionsAction : public ViewerAction
{
public:
    SetPlotOptionsAction(ViewerWindow *);
    virtual ~SetPlotOptionsAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: SetPlotFrameRangeAction
//
// Purpose:
//   This action setes the frames over which a plot exists.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:48:58 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class SetPlotFrameRangeAction : public ViewerAction
{
public:
    SetPlotFrameRangeAction(ViewerWindow *);
    virtual ~SetPlotFrameRangeAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: DeletePlotKeyframeAction
//
// Purpose:
//   Deletes a plot keyframe.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:48:58 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class DeletePlotKeyframeAction : public ViewerAction
{
public:
    DeletePlotKeyframeAction(ViewerWindow *);
    virtual ~DeletePlotKeyframeAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: MovePlotKeyframeAction
//
// Purpose:
//   Moves a plot keyframe.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:48:59 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class MovePlotKeyframeAction : public ViewerAction
{
public:
    MovePlotKeyframeAction(ViewerWindow *);
    virtual ~MovePlotKeyframeAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: SetPlotDatabaseStateAction
//
// Purpose:
//   Sets a database keyframe for a plot.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:49:00 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class SetPlotDatabaseStateAction : public ViewerAction
{
public:
    SetPlotDatabaseStateAction(ViewerWindow *);
    virtual ~SetPlotDatabaseStateAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: DeletePlotDatabaseKeyframeAction
//
// Purpose:
//   Deletes a database keyframe.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:49:02 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class DeletePlotDatabaseKeyframeAction : public ViewerAction
{
public:
    DeletePlotDatabaseKeyframeAction(ViewerWindow *);
    virtual ~DeletePlotDatabaseKeyframeAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: MovePlotDatabaseKeyframeAction
//
// Purpose:
//   Moves a database keyframe.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Apr 11 07:49:03 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class MovePlotDatabaseKeyframeAction : public ViewerAction
{
public:
    MovePlotDatabaseKeyframeAction(ViewerWindow *);
    virtual ~MovePlotDatabaseKeyframeAction();
    virtual void Execute();
};

// ****************************************************************************
// Class: CopyPlotAction
//
// Purpose:
//   This action copies the window's active plots.
//
// Notes:      
//
// Programmer: Ellen Tarwater
// Creation:   Fri Sept 28  2007
//
// Modifications:
//   
// ****************************************************************************

class CopyPlotAction : public ViewerAction
{
public:
    CopyPlotAction(ViewerWindow *);
    virtual ~CopyPlotAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: SetPlotDescriptionAction
//
// Purpose:
//   This action sets the plot description
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Oct 20 13:45:41 PDT 2009
//
// Modifications:
//   
// ****************************************************************************

class SetPlotDescriptionAction : public ViewerAction
{
public:
    SetPlotDescriptionAction(ViewerWindow *);
    virtual ~SetPlotDescriptionAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: MovePlotOrderTowardFirstAction
//
// Purpose:
//   This action moves a plot one slot closer to the plot list start.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Oct 20 13:45:41 PDT 2009
//
// Modifications:
//   
// ****************************************************************************

class MovePlotOrderTowardFirstAction : public ViewerAction
{
public:
    MovePlotOrderTowardFirstAction(ViewerWindow *);
    virtual ~MovePlotOrderTowardFirstAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: MovePlotOrderTowardLastAction
//
// Purpose:
//   This action moves a plot one slot closer to the plot list end.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Oct 20 13:45:41 PDT 2009
//
// Modifications:
//   
// ****************************************************************************

class MovePlotOrderTowardLastAction : public ViewerAction
{
public:
    MovePlotOrderTowardLastAction(ViewerWindow *);
    virtual ~MovePlotOrderTowardLastAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: SetPlotOrderToFirstAction
//
// Purpose:
//   This action moves a plot to the plot list start.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Oct 20 13:45:41 PDT 2009
//
// Modifications:
//   
// ****************************************************************************

class SetPlotOrderToFirstAction : public ViewerAction
{
public:
    SetPlotOrderToFirstAction(ViewerWindow *);
    virtual ~SetPlotOrderToFirstAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

// ****************************************************************************
// Class: SetPlotOrderToLastAction
//
// Purpose:
//   This action moves a plot to the plot list start.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Oct 20 13:45:41 PDT 2009
//
// Modifications:
//   
// ****************************************************************************

class SetPlotOrderToLastAction : public ViewerAction
{
public:
    SetPlotOrderToLastAction(ViewerWindow *);
    virtual ~SetPlotOrderToLastAction();
    virtual void Execute();
    virtual bool Enabled() const;
};

#endif


