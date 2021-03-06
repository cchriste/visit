/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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

#include <visit-config.h>
#include <qapplication.h>
#include <XMLEdit.h>

#include <visitstream.h>
#include <stdlib.h>

#include <QFileDialog>

#include <QString>
#include <QTextStream>
QTextStream cOut(stdout), cErr(stderr);
QString Endl("\n");

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
#ifndef __APPLE__
#include <qwindowsstyle.h>
#endif
#endif

// ****************************************************************************
//  Main Function: XMLEditMain()  
//
//  Purpose:
//    initialize and start the main window
//
//  Programmer:  Jeremy Meredith
//  Creation:    September 26, 2002
//
//  Modifications:
//    Brad Whitlock, Fri Aug 15 15:28:27 PST 2003
//    I prevented the style from being set on MacOS X so it looks like a
//    native Mac application.
//
//    Cyrus Harrison, Thu May 15 16:00:46 PDT 200
//    First pass at porting to Qt 4.4.0
//
// ****************************************************************************

int
XMLEditMain( int argc, char **argv )
{
    QApplication *a = new QApplication(argc, argv);
    XMLEdit *w;
    if (argc > 1)
        w = new XMLEdit(argv[1], NULL);
    else
    {
        w = new XMLEdit("untitled.xml", NULL);
    }

    w->show();

    try
    {
        return a->exec();
    }
    catch (const char *s)
    {
        cerr << "ERROR: " << s << endl;
    }
    catch (const QString &s)
    {
        cerr << "ERROR: " << s.toStdString() << endl;
    }
    return -1;
}

// ****************************************************************************
// Method: main/WinMain
//
// Purpose: 
//   The program entry point function.
//
// Programmer: Brad Whitlock
// Creation:   Wed Nov 23 13:15:31 PST 2011
//
// Modifications:
//   
// ****************************************************************************

#if defined(_WIN32) && defined(VISIT_WINDOWS_APPLICATION)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

int WINAPI
WinMain(HINSTANCE hInstance,     // handle to the current instance
        HINSTANCE hPrevInstance, // handle to the previous instance    
        LPSTR lpCmdLine,         // pointer to the command line
        int nCmdShow             // show state of window
)
{
    return XMLEditMain(__argc, __argv);
}
#else
int
main(int argc, char **argv)
{
    return XMLEditMain(argc, argv);
}
#endif
