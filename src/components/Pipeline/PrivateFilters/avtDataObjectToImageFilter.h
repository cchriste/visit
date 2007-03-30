// ************************************************************************* //
//                        avtDataObjectToImageFilter.h                       //
// ************************************************************************* //

#ifndef AVT_DATA_OBJECT_TO_IMAGE_FILTER_H
#define AVT_DATA_OBJECT_TO_IMAGE_FILTER_H
#include <pipeline_exports.h>


#include <avtFilter.h>
#include <avtImageSource.h>


// ****************************************************************************
//  Class: avtDataObjectToImageFilter
//
//  Purpose:
//      A filter that takes in a data object as input and has an image as
//      output.
//
//  Programmer: Hank Childs
//  Creation:   June 4, 2001
//
// ****************************************************************************

class PIPELINE_API avtDataObjectToImageFilter
    : virtual public avtFilter, virtual public avtImageSource
{
  public:
                       avtDataObjectToImageFilter() {;};
    virtual           ~avtDataObjectToImageFilter() {;};
};


#endif


