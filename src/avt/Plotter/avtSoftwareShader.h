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

// ************************************************************************* //
//                               avtSoftwareShader                           //
// ************************************************************************* //

#ifndef AVT_SOFTWARE_SHADER_H
#define AVT_SOFTWARE_SHADER_H
#include <plotter_exports.h>
#include <avtImage.h>


class  LightAttributes;

struct  avtView3D;


// ****************************************************************************
//  Class: avtSoftwareShader
//
//  Purpose:
//      Does shadows in software.  Also handles depth cueing.
//
//  Programmer: Hank Childs
//  Creation:   October 24, 2004
//
//  Modifications:
//    Jeremy Meredith, Fri Oct 29 16:49:43 PDT 2004
//    Added FindLightView.  Removed the "aspect" argument from AddShadows
//    because we can just calculate it based on the image passed in.
//
//    Jeremy Meredith, Wed Aug 29 13:11:37 EDT 2007
//    Added depth cueing.
//
//    Tom Fogal, Mon Jun 16 11:17:53 EDT 2008
//    Added some const qualifications.
//
//    Jeremy Meredith, Fri Apr 30 14:23:19 EDT 2010
//    Added automatic mode for depth cueing.
//
// ****************************************************************************

class PLOTTER_API avtSoftwareShader
{
  public:
    static bool  GetLightDirection(const LightAttributes &, const avtView3D &,
                                   double *);
    static void  AddShadows(avtImage_p, avtImage_p, const avtView3D &,
                            const avtView3D &, double);
    static void  AddDepthCueing(avtImage_p current_image,
                                const avtView3D &current_view,
                                bool autoExtents,
                                const double startPoint[3],
                                const double endPoint[3],
                                unsigned char cuecolor[3]);
    static avtView3D  FindLightView(avtImage_p, const avtView3D &,
                                    const double*,double);
};

#endif


