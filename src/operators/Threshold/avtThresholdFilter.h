// ************************************************************************* //
//  File: avtThresholdFilter.h
// ************************************************************************* //

#ifndef AVT_Threshold_FILTER_H
#define AVT_Threshold_FILTER_H

#include <avtPluginStreamer.h>
#include <ThresholdAttributes.h>

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
// ****************************************************************************

class avtThresholdFilter : public avtPluginStreamer
{
  public:
                         avtThresholdFilter() {;};
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
    virtual vtkDataSet   *ExecuteData(vtkDataSet *, int, std::string);
    virtual void          RefashionDataObjectInfo(void);
    virtual void          PreExecute(void);
};


#endif
