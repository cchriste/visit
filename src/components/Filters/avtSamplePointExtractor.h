// ************************************************************************* //
//                           avtSamplePointExtractor.h                       //
// ************************************************************************* //

#ifndef AVT_SAMPLE_POINT_EXTRACTOR_H
#define AVT_SAMPLE_POINT_EXTRACTOR_H

#include <pipeline_exports.h>

#include <avtDatasetToSamplePointsFilter.h>
#include <avtVolume.h>


class  vtkHexahedron;
class  vtkPyramid;
class  vtkTetra;
class  vtkVoxel;
class  vtkWedge;

class  avtHexahedronExtractor;
class  avtMassVoxelExtractor;
class  avtPyramidExtractor;
class  avtTetrahedronExtractor;
class  avtWedgeExtractor;

class  avtRayFunction;


// ****************************************************************************
//  Class: avtSamplePointExtractor
//
//  Purpose:
//      This is a component that will take an avtDataset as an input and find
//      all of the sample points from that dataset.
//
//  Programmer: Hank Childs
//  Creation:   December 5, 2000
//
//  Modifications:
//
//    Hank Childs, Sat Jan 27 15:09:34 PST 2001
//    Added support for sending cells when doing parallel volume rendering.
//
//    Kathleen Bonnell, Sat Apr 21, 13:09:27 PDT 2001 
//    Added recursive Execute method to walk down input data tree. 
//
//    Hank Childs, Tue Nov 13 15:51:15 PST 2001
//    Remove boolean argument to Extract<Cell> calls since it is no longer
//    necessary when all of the variables are being extracted.
//
//    Hank Childs, Sun Dec 14 11:07:56 PST 2003
//    Added mass voxel extractor.
//
// ****************************************************************************

class PIPELINE_API avtSamplePointExtractor : public avtDatasetToSamplePointsFilter
{
  public:
                              avtSamplePointExtractor(int, int, int);
    virtual                  ~avtSamplePointExtractor();

    virtual const char       *GetType(void)
                                         { return "avtSamplePointExtractor"; };
    virtual const char       *GetDescription(void)
                                         { return "Extracting sample points";};

    void                      RegisterRayFunction(avtRayFunction *rf)
                                         { rayfoo = rf; };
    void                      SendCellsMode(bool);

  protected:
    int                       width, height, depth;
    int                       currentNode, totalNodes;

    avtHexahedronExtractor   *hexExtractor;
    avtMassVoxelExtractor    *massVoxelExtractor;
    avtPyramidExtractor      *pyramidExtractor;
    avtTetrahedronExtractor  *tetExtractor;
    avtWedgeExtractor        *wedgeExtractor;

    bool                      sendCells;
    avtRayFunction           *rayfoo;

    virtual void              Execute(void);
    virtual void              ExecuteTree(avtDataTree_p);
    void                      SetUpExtractors(void);

    inline void               ExtractHex(vtkHexahedron*,vtkDataSet*,int);
    inline void               ExtractVoxel(vtkVoxel *, vtkDataSet *,int);
    inline void               ExtractTet(vtkTetra *, vtkDataSet *, int);
    inline void               ExtractPyramid(vtkPyramid *, vtkDataSet *, int);
    inline void               ExtractWedge(vtkWedge *, vtkDataSet *, int);
};


#endif


