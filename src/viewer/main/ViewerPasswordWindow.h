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

#ifndef VIEWERPASSWORDWINDOW_H
#define VIEWERPASSWORDWINDOW_H
#include <VisItPasswordWindow.h>
#include <set>

// Forward declarations
class ViewerConnectionProgressDialog;

// ****************************************************************************
//  Class:  ViewerPasswordWindow
//
//  Purpose:
//    Main window for the program.
//
//  Programmer:  Jeremy Meredith
//  Creation:    April 25, 2001
//
//  Modifications:
//    Brad Whitlock, Thu Aug 29 10:40:19 PDT 2002
//    I removed the Windows API stuff since this file is not part of the
//    viewer on Windows. I also added a new static function.
//
//    Brad Whitlock, Thu Aug 29 17:50:25 PST 2002
//    I added a userName argument to the getPassword and authenticate methods.
//
//    Jeremy Meredith, Thu May 24 10:57:02 EDT 2007
//    Added support for checking failed SSH tunneling (port forwards).
//
//    Thomas R. Treadway, Mon Oct  8 13:27:42 PDT 2007
//    Backing out SSH tunneling on Panther (MacOS X 10.3)
//
//    Hank Childs, Sat Nov 10 11:31:34 PST 2007
//    Add a new button for changing the username.
//
//    Kathleen Bonnell, Wed Feb 13 14:05:03 PST 2008 
//    Added static methods to retrieve and reset the value of 
//    'needToChangeUsername'. 
//
//    Brad Whitlock, Tue May 27 13:44:12 PDT 2008
//    Qt 4.
//
//    Kathleen Bonnell, Thu Apr 22 18:06:28 MST 2010 
//    Changed return type of getPassword to std::string.
//
//    Brad Whitlock, Tue Jun 12 11:17:23 PDT 2012
//    I made it a subclass of VisItPasswordWindow.
//
// ****************************************************************************

class ViewerPasswordWindow : public VisItPasswordWindow
{
    Q_OBJECT
public:
    ViewerPasswordWindow(QWidget *parent=NULL);
    virtual ~ViewerPasswordWindow();

    // Callback function for RemoteProcess' authentication callback.
    static void authenticate(const char *username, const char* password, const char *host, int fd);

    static void SetConnectionProgressDialog(ViewerConnectionProgressDialog *d);

    static std::set<int> GetFailedPortForwards();
private:
    std::string password(const char *username, const char *host,
                         bool passphrase, VisItPasswordWindow::ReturnCode &ret);

    static ViewerPasswordWindow           *instance;
    static ViewerConnectionProgressDialog *dialog;
    static std::set<int>                   failedPortForwards;
};

#endif
