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

#include "Explorer.h"
#include <SiloView.h>

#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qfiledialog.h>
#include <qapplication.h>

// ****************************************************************************
//  Constructor: Explorer::Explorer
//
//  Programmer:  Jeremy Meredith
//  Creation:    November 12, 2001
//
//  Modifications:
//    Mark C. Miller, Thu Jul 20 15:45:55 PDT 2006
//    Made it deal with failure to construct SiloView
//  
// ****************************************************************************
Explorer::Explorer(const QString &file, QWidget *p, const QString &n)
    : QMainWindow(p,n)
{
    view = new SiloView(file,
                        this, "SiloView");
    if (!view->HasSiloFile())
    {
        delete view;
        view = 0;
    }
    else
    {
        setCentralWidget(view);
        setCaption("Explorer: "+file);

        QPopupMenu *filemenu = new QPopupMenu( this );
        menuBar()->insertItem(tr("&File"),filemenu);
        filemenu->insertItem( "&Open",  this, SLOT(open()),  CTRL+Key_O );
        filemenu->insertSeparator();
        filemenu->insertItem( "E&xit", this, SLOT(close()),  CTRL+Key_X );
    }
}

// ****************************************************************************
//  Destructor:  Explorer::~Explorer
//
//  Programmer:  Jeremy Meredith
//  Creation:    May 17, 2004
//
//  Modifications:
//    Mark C. Miller,Thu Jul 20 15:45:55 PDT 2006
//    Added deletion of view
// ****************************************************************************
Explorer::~Explorer()
{
    if (view) delete view;
    view = 0;
}

// ****************************************************************************
//  Method:  Explorer::open
//
//  Purpose:
//    Open a new file.
//
//  Programmer:  Jeremy Meredith
//  Creation:    November 12, 2001
//
//  Modifications:
//    Jeremy Meredith, Mon May 17 11:50:08 PDT 2004
//    Change the window caption when opening a new file.
//
// ****************************************************************************
void
Explorer::open()
{
    QString file = 
        QFileDialog::getOpenFileName(QString(),
                                     "Silo files (*.silo *.root *.pdb);;"
                                     "All files (*)",
                                     NULL, "SiloOpen", "Open file...");
    if (file.isNull())
        return;

    setCaption("Explorer: "+file);
    view->Set(file);
}
