// ************************************************************************* //
//                                 avtHistogramPlot.h                             //
// ************************************************************************* //

#ifndef AVT_Histogram_PLOT_H
#define AVT_Histogram_PLOT_H


#include <avtLegend.h>
#include <avtPlot.h>
#include <avtSurfaceAndWireframeRenderer.h>

#include <HistogramAttributes.h>

class     avtHistogramFilter;
class     avtUserDefinedMapper;
class     vtkProperty;


// ****************************************************************************
//  Class:  avtHistogramPlot
//
//  Purpose:
//      A concrete type of avtPlot, this is the Histogram plot.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Thu Jun 26 10:33:56 PDT 2003
//
// ****************************************************************************

class avtHistogramPlot : public avtLineDataPlot
{
  public:
                                avtHistogramPlot();
    virtual                    ~avtHistogramPlot();

    virtual const char         *GetName(void) { return "HistogramPlot"; };

    static avtPlot             *Create();
    virtual void                ReleaseData(void);

    virtual void                SetAtts(const AttributeGroup*);
    virtual bool                SetForegroundColor(const double *);

  protected:
    HistogramAttributes              atts;

    avtSurfaceAndWireframeRenderer_p renderer;
    avtUserDefinedMapper            *mapper;
    avtHistogramFilter              *HistogramFilter;
    avtFilter                       *amountFilter;
    vtkProperty                     *property; 
    double                           fgColor[3];

    virtual avtMapper          *GetMapper(void);
    virtual avtDataObject_p     ApplyOperators(avtDataObject_p);
    virtual avtDataObject_p     ApplyRenderingTransformation(avtDataObject_p);
    virtual void                CustomizeBehavior(void);

    virtual avtLegend_p         GetLegend(void) { return NULL; };
};


#endif
