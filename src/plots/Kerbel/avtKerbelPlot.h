// ************************************************************************* //
//                                 avtKerbelPlot.h                             //
// ************************************************************************* //

#ifndef AVT_Kerbel_PLOT_H
#define AVT_Kerbel_PLOT_H


#include <avtLegend.h>
#include <avtPlot.h>

#include <KerbelAttributes.h>
#include <avtVariableMapper.h>

class     avtVariableLegend;
class     avtKerbelFilter;
class     avtLookupTable;


// ****************************************************************************
//  Class:  avtKerbelPlot
//
//  Purpose:
//      A concrete type of avtPlot, this is the Kerbel plot.
//
//  Programmer: ahern -- generated by xml2info
//  Creation:   Tue Dec 18 16:23:56 PST 2001
//
//  Modifications:
//    Kathleen Bonnell, Thu Sep 18 13:44:50 PDT 2003
//    Moved SetOpacity to C file.
//
// ****************************************************************************

class avtKerbelPlot : public avtSurfaceDataPlot
{
  public:
                                avtKerbelPlot();
    virtual                    ~avtKerbelPlot();

    virtual const char         *GetName(void) { return "KerbelPlot"; };

    static avtPlot             *Create() { return new avtKerbelPlot; };

    void                        SetOpacity(float o);
    virtual bool                SetColorTable(const char *);
    void                        SetLegend(bool);
    void                        SetVarName(const char *);
    virtual void                SetAtts(const AttributeGroup*);

  protected:
    KerbelAttributes              atts;

    bool                          colorsInitialized;
    avtLookupTable               *avtLUT;

    avtVariableMapper            *varMapper;
    avtVariableLegend            *varLegend;
    avtLegend_p                  varLegendRefPtr;
    avtKerbelFilter              *filter;

    virtual avtMapper          *GetMapper(void) { return varMapper; };
    virtual avtLegend_p         GetLegend(void) { return varLegendRefPtr; };
    virtual avtDataObject_p     ApplyRenderingTransformation(avtDataObject_p p) {return p;};
    virtual avtDataObject_p     ApplyOperators(avtDataObject_p);
    virtual void                CustomizeBehavior(void);
    virtual void                CustomizeMapper(avtDataObjectInformation &);
};


#endif
