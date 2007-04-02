// ************************************************************************* //
//  File: avtThresholdFilter.h
// ************************************************************************* //

#ifndef AVT_Threshold_FILTER_H
#define AVT_Threshold_FILTER_H

#include <avtPluginStructuredChunkStreamer.h>
#include <ThresholdAttributes.h>

#include <avtGhostData.h>

class     vtkDataSet;


// ****************************************************************************
//  Class: avtThresholdFilter
//
//  Purpose:
//      A plugin operator for Threshold.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Tue Oct 23 16:38:18 PST 2001
//
//  Modifications:
//
//    Hank Childs, Sat Mar 19 10:18:52 PST 2005
//    Add support for structured chunking.
//
//    Hank Childs, Sun Mar 27 11:34:04 PST 2005
//    Inherit from new base type that supports structured chunking.
//
//    Hank Childs, Tue Sep 13 09:07:05 PDT 2005
//    Add support for PointsOnly mode.
//
// ****************************************************************************

class avtThresholdFilter : public avtPluginStructuredChunkStreamer
{
  public:
                         avtThresholdFilter();
    virtual             ~avtThresholdFilter() {;};

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtThresholdFilter"; };
    virtual const char  *GetDescription(void)
                             { return "Thresholding"; };

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    ThresholdAttributes   atts;

    virtual avtPipelineSpecification_p
                          PerformRestriction(avtPipelineSpecification_p);

    virtual vtkDataSet   *ProcessOneChunk(vtkDataSet *, int, std::string,bool);
    virtual void          GetAssignments(vtkDataSet *, const int *,
                      std::vector<avtStructuredMeshChunker::ZoneDesignation>&);

    virtual void          RefashionDataObjectInfo(void);
    virtual void          PreExecute(void);

    vtkDataArray         *GetThresholdVariable(vtkDataSet *, bool &);
    vtkDataSet           *ThresholdToPointMesh(vtkDataSet *);
};


#endif
