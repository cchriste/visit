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

#ifndef CELL_RECONSTRUCTOR_H
#define CELL_RECONSTRUCTOR_H

#include <ZooMIR.h>
#include <VisItArray.h>

// ****************************************************************************
//  Class:  CellReconstructor
//
//  Purpose:
//    Reconstructs a cell from the given mesh/material.
//
//  Programmer:  Jeremy Meredith
//  Creation:    September 15, 2003
//
//  Modifications:
//    Jeremy Meredith, Thu Aug 18 18:02:38 PDT 2005
//    I was able to re-use most of this class for a new isovolume based
//    reconstruction algorithm.  Everything stayed except I made
//    ReconstructCell pure-virtual, and I needed to keep track of whether
//    or not edge points were shared across materials.
//
//    Mark C. Miller, Thu Feb  9 21:06:10 PST 2006
//    Renamed Array class to VisItArray to avoid name collisions with
//    third-party libs
//
//    Jeremy Meredith, Fri Feb 13 10:56:43 EST 2009
//    Allowed ReconstructCell to take a new argument where it will, if
//    desired, output the volume fractions for the reconstructed materials.
//    Also, added helper function to calculate volume or area.
//
//    Jeremy Meredith, Tue Jun 18 11:56:22 EDT 2013
//    Output actual volumes/areas, not VF's, and return total vol/area, 
//    in ReconstructCell.
//
// ****************************************************************************
class CellReconstructor
{
  public:
    CellReconstructor(vtkDataSet*, avtMaterial*, ResampledMat&, int, int, bool,
                      MIRConnectivity&, ZooMIR&);
    virtual ~CellReconstructor();
    virtual double ReconstructCell(int, int, int, vtkIdType*, double*) = 0;

  protected:
    vtkDataSet                             *mesh;
    avtMaterial                            *mat;
    avtMaterial                            *origMat;
    ResampledMat                           &rm;
    int                                     nPoints;
    int                                     nCells;
    MIRConnectivity                        &conn;
    ZooMIR                                 &mir;
    int                                     nMaterials;

    static double CalculateVolumeOrAreaHelper(int celltype,double coords[][3]);

  protected:
    int           cellid;
    int           celltype;
    vtkIdType    *ids;
    int           nids; 

    int           nodeIndices[MAX_NODES_PER_ZONE];
    int           nodeIndexLimit[MAX_NODES_PER_ZONE];
    float        *nodeVFs[MAX_NODES_PER_ZONE];
    float         newVFs[MAX_NODES_PER_ZONE];
    float         vfDiff[MAX_NODES_PER_ZONE];
    int          *mix_index;
    int           interpIDs[4];
    float         interpVFs[4];

    bool            allMaterialsSharePoints;
    int             startIndex;
    unsigned char  *splitCase;
    int             numOutput;
    typedef int     edgeIndices[2];
    edgeIndices    *vertices_from_edges;

    struct TempCell
    {
        int mat;
        int celltype;
        int nnodes;
        int mix_index;
        int ids[MAX_NODES_PER_ZONE];
        int localids[MAX_NODES_PER_ZONE];
        float vfs[MAX_NODES_PER_ZONE];
    };

    ZooMIR::EdgeHashTable  edges;
    VisItArray<TempCell>   outlist;
    VisItArray<TempCell>   tmplist;

    void CreateCentroidPoint(TempCell&, int, int, int, int, int, int*);
    void CreateOutputShape(TempCell&, TempCell&, int, int, int*, int);

};

#endif
