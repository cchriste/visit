/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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
//                        vtkVolumeFromCSGVolume.h                           //
// ************************************************************************* //

#ifndef VTK_VOLUME_FROM_CSG_VOLUME_H
#define VTK_VOLUME_FROM_CSG_VOLUME_H

#include <visit_vtk_exports.h>
#include <vtkVolumeFromVolume.h>
#include <FixedLengthBitField.h>

#include <vector>


class vtkCellData;
class vtkPointData;
class vtkPolyData;
class vtkUnstructuredGrid;


// ****************************************************************************
//  Class: vtkVolumeFromCSGVolume
//
//  Purpose:
//      This is a special type of vtkVolumeFromVolume, where each added
//      cell belongs to one half of a binary partition of the original
//      data set.  You specify which half when you add the shape, and
//      when constructing the data set, tell it the index (into a bit array)
//      corresponding to the current binary partition function.
//      E.g. on the first pass, you specify bit #0 and "in" cells get the tag
//      0x0 and "out" cels get the tag 0x1.  On the second pass, you specify
//      bit #1, and "in" cells get tagged either 0x0 or 0x1, and "out" cells
//      get tagged either 0x2 or 0x3.  And so on, up to the width of the
//      bit field used here.
//
//  Programmer: Eric Brugger
//  Creation:   July 23, 2012
//
//  Modifications:
//
// ****************************************************************************

class VISIT_VTK_API vtkVolumeFromCSGVolume : public vtkVolumeFromVolume
{
  protected:


  public:
                      vtkVolumeFromCSGVolume(int nPts,
                                             int ptSizeGuess);
    virtual          ~vtkVolumeFromCSGVolume() { ; };

    void              ConstructDataSet(vtkCellData *,
                                       vtkUnstructuredGrid *, float *, int,
                                       std::vector<FixedLengthBitField<64> > *newTags);

    using vtkVolumeFromVolume::AddCentroidPoint;
    using vtkDataSetFromVolume::AddPoint;

    void           AddHex(int z, int v0, int v1, int v2, int v3,
                          int v4, int v5, int v6, int v7,
                          FixedLengthBitField<64> &tag)
    {
        vtkVolumeFromVolume::AddHex(z, v0, v1, v2, v3, v4, v5, v6, v7);
        hexTags.push_back(tag);
    }
        
    void           AddWedge(int z,int v0,int v1,int v2,int v3,int v4,int v5,
                            FixedLengthBitField<64> &tag)
    {
        vtkVolumeFromVolume::AddWedge(z, v0, v1, v2, v3, v4, v5);
        wedgeTags.push_back(tag);
    }
    void           AddPyramid(int z, int v0, int v1, int v2, int v3, int v4,
                              FixedLengthBitField<64> &tag)
    {
        vtkVolumeFromVolume::AddPyramid(z, v0, v1, v2, v3, v4); 
        pyramidTags.push_back(tag);
    }
    void           AddTet(int z, int v0, int v1, int v2, int v3,
                          FixedLengthBitField<64> &tag)
    {
        vtkVolumeFromVolume::AddTet(z, v0, v1, v2, v3);
        tetTags.push_back(tag);
    }
    void           AddQuad(int z, int v0, int v1, int v2, int v3,
                           FixedLengthBitField<64> &tag)
    {
        vtkVolumeFromVolume::AddQuad(z, v0, v1, v2, v3);
        quadTags.push_back(tag);
    }
    void           AddTri(int z, int v0, int v1, int v2,
                          FixedLengthBitField<64> &tag)
    {
        vtkVolumeFromVolume::AddTri(z, v0, v1, v2);
        triTags.push_back(tag);
    }
    void           AddLine(int z, int v0, int v1,
                           FixedLengthBitField<64> &tag)
    {
        vtkVolumeFromVolume::AddLine(z, v0, v1);
        lineTags.push_back(tag);
    }
    void           AddVertex(int z, int v0,
                             FixedLengthBitField<64> &tag)
    {
        vtkVolumeFromVolume::AddVertex(z, v0);
        vertexTags.push_back(tag);
    }


    void                    InitTraversal();
    int                     GetNumberOfCells() const;
    void                    NextCell();
    const vtkIdType        *GetCell();
    FixedLengthBitField<64> GetTag();
    int                     GetCellSize() const;
    int                     GetCellVTKType() const;
    void                    SetId(int);
    void                    SetTagBit(int);
    void                    InvalidateCell();

    void                    UpdatePoints(std::vector<float> &);

  protected:
    std::vector<FixedLengthBitField<64> >   hexTags;
    std::vector<FixedLengthBitField<64> >   wedgeTags;
    std::vector<FixedLengthBitField<64> >   pyramidTags;
    std::vector<FixedLengthBitField<64> >   tetTags;
    std::vector<FixedLengthBitField<64> >   quadTags;
    std::vector<FixedLengthBitField<64> >   triTags;
    std::vector<FixedLengthBitField<64> >   lineTags;
    std::vector<FixedLengthBitField<64> >   vertexTags;
    std::vector<FixedLengthBitField<64> >  *shapeTags[8];
    std::vector<FixedLengthBitField<64> >  *curTags;

    int                     ishape;
    int                     shapeCnt[8];
    ShapeList              *curShapes;
    int                     curShapeCnt;
    int                     curShapeSize;
    int                     curShapeVTKType;
    int                     curShape;
};

#endif
