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

#include <visit_writer.h>

int
main(int argc, char *argv[])
{
#define NX 4
#define NY 5
#define NZ 3
    /* Rectilinear mesh coordinates. */
    float x[] = {0., 1., 2.5, 5.};
    float y[] = {0., 2., 2.25, 2.55,  5.};
    float z[] = {0., 1., 3.};
    int dims[] = {NX, NY, NZ};
    int ndims = 3;
    /* Zonal and Nodal variable data. */
    float zonal[NZ-1][NY-1][NX-1], nodal[NZ][NY][NX];
    float zonalvec[NZ-1][NY-1][NX-1][3], nodalvec[NZ][NY][NX][3];
    /* Info about the variables to pass to visit_writer. */
    int nvars = 4;
    int vardims[] = {1, 1, 3, 3};
    int centering[] = {0, 1, 0, 1};
    const char *varnames[] = {"zonal", "nodal", "zonalvec", "nodalvec"};
    float *vars[] = {(float*)zonal, (float*)nodal, (float*)zonalvec, (float*)nodalvec};

    /* Create 2 zonal variables; 1 scalar, 1 vector. */
    int i,j,k,index = 0;
    for(k = 0; k < NZ-1; ++k)
        for(j = 0; j < NY-1; ++j)
            for(i = 0; i < NX-1; ++i, ++index)
            {
                zonal[k][j][i] = (float)index;

                zonalvec[k][j][i][0] = 1.f;
                zonalvec[k][j][i][1] = 0.f;
                zonalvec[k][j][i][2] = 0.f;
            }

    /* Create 2 nodal variables; 1 scalar, 1 vector. */
    index = 0;
    for(k = 0; k < NZ; ++k)
        for(j = 0; j < NY; ++j)
            for(i = 0; i < NX; ++i, ++index)
            {
                nodal[k][j][i] = (float)index;

                nodalvec[k][j][i][0] = 0.f;
                nodalvec[k][j][i][1] = 1.f;
                nodalvec[k][j][i][2] = 0.f;
            }

    /* Pass the data to visit_writer to write a VTK file.*/
    write_rectilinear_mesh("vwrect3d.vtk", 0, dims, x, y, z, nvars,
        vardims, centering, varnames, vars);

    return 0;
}
