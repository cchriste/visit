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

#ifndef AVT_OPENGL_ATOM_TEXTURER_3D_H
#define AVT_OPENGL_ATOM_TEXTURER_3D_H

// ****************************************************************************
// Class: avtOpenGLAtomTexturer3D
//
// Purpose:
//   Performs sphere shading on imposter quads using different methods.
//
// Notes:      There are 2 methods: texture-based, shader-based.
//
// Programmer: Brad Whitlock
// Creation:   Tue Mar 28 09:56:02 PDT 2006
//
// Modifications:
//   Brad Whitlock, Fri Apr 7 11:26:48 PDT 2006
//   Added SetHint method.
//
//   John Schreiner, Fri Feb 12 19:19:34 MST 2010
//   Removed width/height hints that aren't required anymore.
//
//   Jeremy Meredith, Thu Apr 22 14:36:12 EDT 2010
//   Renamed everything to have a "3D" in the name.
//
// ****************************************************************************

class avtOpenGLAtomTexturer3D
{
public:
    avtOpenGLAtomTexturer3D();
    virtual ~avtOpenGLAtomTexturer3D();

    void BeginSphereTexturing();
    void EndSphereTexturing();

    static const int HINT_SET_DEPTH;

    void SetHint(int hint, int value);

    typedef enum {TextureMode, ShaderMode} TexMode_t;
    TexMode_t GetMode();

private:
    bool modeDetermined;
    TexMode_t mode;
    void *tData;
    void *sData;
};


#endif
