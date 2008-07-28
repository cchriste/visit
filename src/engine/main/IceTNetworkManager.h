/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
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

#ifndef ICET_NETWORK_MANAGER_H
#define ICET_NETWORK_MANAGER_H

#include <NetworkManager.h>
#include <GL/ice-t.h>

// ****************************************************************************
//  Class: IceTNetworkManager
//
//  Purpose:
//      NetworkManager which uses IceT for rendering/readback of image data.
//
//  Programmer: Tom Fogal
//  Creation:   June 17, 2008
//
//  Modifications:
//
//    Tom Fogal, Tue Jun 24 13:27:48 EDT 2008
//    Defined `Readback' function.
//
//    Tom Fogal, Mon Jul 14 12:27:23 PDT 2008
//    Override parent's timer information.
//
//    Tom Fogal, Wed Jul 16 12:59:37 EDT 2008
//    Oops, destructor should be virtual.
//
//    Tom Fogal, Sat Jul 26 23:07:15 EDT 2008
//    Override RenderGeometry for a potential IceT-only optimization.
//
// ****************************************************************************

class IceTNetworkManager: public NetworkManager
{
 public:
               IceTNetworkManager(void);
    virtual   ~IceTNetworkManager(void);

    void       TileLayout(size_t width, size_t height) const;

    virtual avtDataObjectWriter_p
               Render(intVector networkIds, bool getZBuffer, int annotMode,
                      int windowID, bool leftEye);
    void       RealRender(); /// OpenGL calls sourced from here

 protected:

    virtual avtImage_p RenderGeometry();
    virtual avtImage_p Readback(const VisWindow * const, bool) const;
    virtual void       StopTimer(int windowID);

 private:

    void  VerifyColorFormat() const;

 private:
    IceTCommunicator comm;
    IceTContext context;
};

#endif /* ICET_NETWORK_MANAGER_H */
