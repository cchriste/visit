/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
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

#ifndef QUATERNION_H
#define QUATERNION_H
#include <math_exports.h>

#include <avtVector.h>
#include <avtMatrix.h>

// ****************************************************************************
//  Purpose:
//    Build a quaternion using a virtual trackball. The virtual trackball
//    is based on a Witch of Agnesi, which is a sphere unioned with a
//    plane.  The size of the sphere should really be based on the distance
//    from the center of rotation to the point on the object underneath the
//    mouse.  That point would track the mouse as closely as possible.
//    The witch of agnesi is parameterized with these 3 parameters:
//
//        AGNESI_RADIUS  0.8    The Z value at x = y = 0.0
//        COMPRESSION    3.5    Multipliers for x and y
//        SPHERE_RADIUS  1.0    Radius of the sphere embedded in the plane
//
// ****************************************************************************

class MATH_API avtQuaternion
{
  public:
    avtQuaternion();
    avtQuaternion(const avtVector&, double);
    avtQuaternion(const avtQuaternion&);
    void operator=(const avtQuaternion&);

    avtMatrix CreateRotationMatrix();

    double norm();
    void normalize();
  private:
    double x;
    double y;
    double z;
    double s;
};

#endif

