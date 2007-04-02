// ************************************************************************* //
//  File: avtSmoothFilter.h
// ************************************************************************* //

#ifndef AVT_Smooth_FILTER_H
#define AVT_Smooth_FILTER_H


#include <avtPluginStreamer.h>
#include <SmoothOperatorAttributes.h>


class vtkDataSet;


// ****************************************************************************
//  Class: avtSmoothFilter
//
//  Purpose:
//      A plugin operator for Smooth.
//
//  Programmer: childs -- generated by xml2avt
//  Creation:   Sun Aug 14 11:59:58 PDT 2005
//
// ****************************************************************************

class avtSmoothFilter : public avtPluginStreamer
{
  public:
                         avtSmoothFilter();
    virtual             ~avtSmoothFilter();

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtSmoothFilter"; };
    virtual const char  *GetDescription(void)
                             { return "Smooth"; };

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    SmoothOperatorAttributes   atts;
    bool                       issuedWarning;

    virtual vtkDataSet   *ExecuteData(vtkDataSet *, int, std::string);
    virtual void          PreExecute(void);
};


#endif
