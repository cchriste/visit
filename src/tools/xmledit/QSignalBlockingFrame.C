/*****************************************************************************
*
* Copyright (c) 2000 - 2012, Lawrence Livermore National Security, LLC
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

#include <QSignalBlockingFrame.h>

#include <vector>
using std::vector;
// ****************************************************************************
//  Constructor:  QSignalBlockingFrame::QSignalBlockingFrame
//
//  Programmer:  Kevin Bensema
//  Creation:    September 11, 2012
//
// ****************************************************************************
QSignalBlockingFrame::QSignalBlockingFrame(QWidget* parent)
  : QFrame(parent) 
{
  ;
}

// ****************************************************************************
//  Method:  QSignalBlockingFrame::BlockAllSignals
//
//  Purpose:
//    Blocks/unblocks signals to the widgets.  This lets them get
//    updated by changes in state without affecting the state.
//    No modifications are needed to this function when new widgets are added.
//    Works via depth-first search, and calls blockSignals on all widgets.
//
//  Arguments:
//    block      whether to block (true) or unblock (false) signals
//
//  Programmer:  Kevin Bensema
//  Creation:    September 11, 2012
//
//  Modifications:
//
// ****************************************************************************
void
QSignalBlockingFrame::BlockAllSignals(bool block)
{
  vector<QObjectList> recursionStack;
  recursionStack.push_back(children());
  while(recursionStack.empty() == false)
  {
    QObjectList listOfChildren = recursionStack.back();
    recursionStack.pop_back();

    // iterate over the list of children, calling blockSignals on each.
    // if a child has children of its own, add them to the stack.
    for(QObjectList::Iterator iter = listOfChildren.begin();
        iter != listOfChildren.end(); ++iter)
    {
      QObject *object = *iter;
      object->blockSignals(block);
      QObjectList newChildren = object->children();
      if(newChildren.empty() == false)
      {
        recursionStack.push_back(newChildren);
      }
    }
  }

  return;
}

