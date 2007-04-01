// ************************************************************************* //
//  File: avtIsovolumeFilter.h
// ************************************************************************* //

#ifndef AVT_Isovolume_FILTER_H
#define AVT_Isovolume_FILTER_H


#include <avtPluginStreamer.h>
#include <IsovolumeAttributes.h>


class vtkDataSet;

// ****************************************************************************
//  Class: avtIsovolumeFilter
//
//  Purpose:
//      A plugin operator for Isovolume.
//
//  Programmer: meredith -- generated by xml2info
//  Creation:   Fri Jan 30 14:50:21 PST 2004
//
//  Modifications:
//    Jeremy Meredith, Mon Feb 16 19:12:11 PST 2004
//    Added RefashionDataObjectInfo.  This was needed for correct support
//    on various mesh types.
//
//    Jeremy Meredith, Thu May  6 11:37:47 PDT 2004
//    Split some code from ExecuteData into a new function to avoid
//    code duplication.
//
// ****************************************************************************

class avtIsovolumeFilter : public avtPluginStreamer
{
  public:
                         avtIsovolumeFilter();
    virtual             ~avtIsovolumeFilter();

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtIsovolumeFilter"; };
    virtual const char  *GetDescription(void)
                             { return "Isovolume"; };

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    IsovolumeAttributes   atts;

    virtual vtkDataSet   *ExecuteData(vtkDataSet *, int, std::string);
    virtual void          RefashionDataObjectInfo(void);

  private:
    virtual vtkDataSet   *ExecuteSingleClip(vtkDataSet *, float, bool);
};


#endif
