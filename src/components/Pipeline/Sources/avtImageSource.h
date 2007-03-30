// ************************************************************************* //
//                            avtImageSource.h                               //
// ************************************************************************* //

#ifndef AVT_IMAGE_SOURCE_H
#define AVT_IMAGE_SOURCE_H
#include <pipeline_exports.h>


#include <avtDataObjectSource.h>
#include <avtImage.h>


class     vtkImageData;

class     avtImageRepresentation;


// ****************************************************************************
//  Class: avtImageSource
//
//  Purpose:
//      A data object source whose data object is an image.
//
//  Programmer: Hank Childs
//  Creation:   June 1, 2001
//
// ****************************************************************************

class PIPELINE_API avtImageSource : virtual public avtDataObjectSource
{
  public:
                                avtImageSource();
    virtual                    ~avtImageSource() {;};

    virtual avtDataObject_p     GetOutput(void);
    vtkImageData               *GetVTKOutput(void);

    avtImage_p                  GetTypedOutput(void) { return image; };

  protected:
    avtImage_p                  image;

    void                        SetOutput(avtImage_p);

    avtImageRepresentation     &GetOutputImage(void);
    void                        SetOutputImage(const avtImageRepresentation &);
};


#endif


