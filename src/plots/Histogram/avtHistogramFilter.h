// ************************************************************************* //
//                           avtHistogramFilter.h                            //
// ************************************************************************* //

#ifndef AVT_Histogram_FILTER_H
#define AVT_Histogram_FILTER_H


#include <avtStreamer.h>

#include <HistogramAttributes.h>


// ****************************************************************************
//  Class: avtHistogramFilter
//
//  Purpose:
//      This operator is the implied operator associated with an Histogram plot.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Thu Jun 26 09:04:54 PDT 2003
//
//  Modifications:
//
//    Hank Childs, Sat Oct 18 12:14:32 PDT 2003
//    Added PerformRestriction to disable ghost zone generation.
//
// ****************************************************************************

class avtHistogramFilter : public avtStreamer
{
  public:
                              avtHistogramFilter();
    virtual                  ~avtHistogramFilter();

    virtual const char       *GetType(void)   { return "avtHistogramFilter"; };
    virtual const char       *GetDescription(void)
                                  { return "Performing Histogram"; };

    void                      SetAttributes(const HistogramAttributes &);
   
  protected:
    HistogramAttributes       atts;
    float                    *bins;
    double                    workingMin;
    double                    workingMax;
    int                       workingNumBins;

    virtual vtkDataSet       *ExecuteData(vtkDataSet *, int, std::string);
    virtual void              RefashionDataObjectInfo(void);
    virtual void              PreExecute();
    virtual void              PostExecute();
    virtual avtPipelineSpecification_p
                              PerformRestriction(avtPipelineSpecification_p);
};


#endif


