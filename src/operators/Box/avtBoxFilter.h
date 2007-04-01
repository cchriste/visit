// ************************************************************************* //
//  File: avtBoxFilter.h
// ************************************************************************* //

#ifndef AVT_Box_FILTER_H
#define AVT_Box_FILTER_H

#include <avtPluginStreamer.h>
#include <BoxAttributes.h>

class     vtkDataSet;
class     vtkRectilinearGrid;
class     vtkUnstructuredGrid;


// ****************************************************************************
//  Class: avtBoxFilter
//
//  Purpose:
//      A plugin operator for Box.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Mon Nov 12 16:57:31 PST 2001
//
//  Modifications:
//
//    Mark C. Miller, Tue Sep 28 19:57:42 PDT 2004
//    Added selection id and PerformRestriction implementation
//
//    Hank Childs, Sun Apr 24 11:11:46 PDT 2005
//    Add special support for rectilinear grids. ['6155]
//
// ****************************************************************************

class avtBoxFilter : public avtPluginStreamer
{
  public:
                         avtBoxFilter();
    virtual             ~avtBoxFilter();

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtBoxFilter"; };
    virtual const char  *GetDescription(void)
                             { return "Finding cells within a box"; };

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    BoxAttributes   atts;
    int             selID;

    virtual vtkDataSet   *ExecuteData(vtkDataSet *, int, std::string);
    vtkRectilinearGrid   *RectilinearExecute(vtkRectilinearGrid *);
    vtkUnstructuredGrid  *GeneralExecute(vtkDataSet *);
    virtual void          RefashionDataObjectInfo(void);
    virtual avtPipelineSpecification_p
                            PerformRestriction(avtPipelineSpecification_p);

};


#endif
