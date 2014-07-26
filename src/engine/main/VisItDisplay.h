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

#ifndef VISIT_VISIT_DISPLAY_H
#define VISIT_VISIT_DISPLAY_H

#include <string>
#include <vector>
#include <engine_main_exports.h>

// ****************************************************************************
//  Class:  VisItDisplay
//
//  Purpose:
//    Interface for initializing a display to utilize in rendering.
//
//  Programmer:  Tom Fogal
//  Creation:    September 1, 2008
//
//  Modifications:
//
//    Tom Fogal, Mon Sep  1 15:11:06 EDT 2008
//    Add a method to create the appropriate display type (factory).
//
//    Tom Fogal, Tue May 25 16:08:39 MDT 2010
//    Made ::Connect return a bool, to detect errors.
//
//    Brad Whitlock, Mon Oct 10 11:37:25 PDT 2011
//    Added GetDisplayType method.
//
// ****************************************************************************

class ENGINE_MAIN_API VisItDisplay
{
  public:
    virtual ~VisItDisplay() {}

    typedef enum {
        D_MESA, // mesa based SW rendering
        D_X     // utilize an X server for HW rendering
    } DisplayType;

    /// Initializes a VisItDisplay.  Rendering is only possible after
    /// `Initialize' and before `Teardown'.  `Initialize' must be called prior
    /// to any rendering commands.  Multiple calls to `Initialize' without
    /// intervening `Teardown's result in undefined behavior.
    /// Initialization may fail; success is given via the return value.
    ///
    /// display the display string to utilize.
    /// args    user arguments appropriate for configuring the VisIt display.
    /// returns: Display initialization status.
    virtual bool   Initialize(std::string display,
                              const std::vector<std::string> &args) = 0;

    /// Associates this process with the previously-`Initialize'd VisIt Display.
    /// Undefined if this VisItDisplay has not been initialized.
    virtual bool   Connect() = 0;

    /// Closes down a VisItDisplay.  You may assume that this will be called by
    /// the destructor.  It is safe to call Teardown more than once.  It is not
    /// safe to use `Teardown' before `Initialize'.  This implicitly breaks any
    /// previous connections.
    virtual void   Teardown() = 0;

    // return the display type.
    virtual DisplayType GetDisplayType() const = 0;
};

// Methods helpful in argument vector creation.
// If these were to be generalized, they would more appropriately end up in
// `StringHelpers'
std::string display_format(std::string s, size_t node, size_t display);
std::vector<std::string> split(std::string, size_t, size_t);

namespace VDisplay {
    /// Creates a display type based on the given parameter.  It is the
    /// caller's responsibility to deallocate the display via `delete'.
    VisItDisplay *Create(VisItDisplay::DisplayType t);
};

#endif /* VISIT_VISIT_DISPLAY_H */
