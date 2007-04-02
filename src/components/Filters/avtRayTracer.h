// ************************************************************************* //
//                                 avtRayTracer.h                            //
// ************************************************************************* //

#ifndef AVT_RAY_TRACER_H
#define AVT_RAY_TRACER_H

#include <filters_exports.h>

#include <avtDatasetToImageFilter.h>
#include <avtViewInfo.h>

class   avtRayFunction;


// ****************************************************************************
//  Class: avtRayTracer
//
//  Purpose:
//      Performs ray tracing, taking in a dataset as a source and has an
//      image as an output.
//
//  Programmer: Hank Childs
//  Creation:   November 27, 2000
//
//  Modifications:
//
//    Hank Childs, Mon Jan  8 16:52:26 PST 2001
//    Added "Get" functions.
//
//    Hank Childs, Sat Feb  3 20:37:01 PST 2001
//    Removed pixelizer and added mechanism to change background color.
//
//    Hank Childs, Tue Feb 13 15:15:50 PST 2001
//    Added ability to insert an opaque image into the rendering.
//
//    Brad Whitlock, Wed Dec 5 11:13:18 PDT 2001
//    Added gradient backgrounds.
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined destructor definition to .C file because certain compilers
//    have problems with them.
//
//    Hank Childs, Sun Dec  4 18:00:55 PST 2005
//    Add method that estimates number of stages.
//
//    Hank Childs, Mon Jan 16 11:11:47 PST 2006
//    Add support for kernel based sampling.
//
// ****************************************************************************

class AVTFILTERS_API avtRayTracer : public avtDatasetToImageFilter
{
  public:
                          avtRayTracer();
    virtual              ~avtRayTracer();

    virtual const char   *GetType(void) { return "avtRayTracer"; };
    virtual const char   *GetDescription(void) { return "Ray tracing"; };
    virtual void          ReleaseData(void);

    void                  SetView(const avtViewInfo &);

    static int            GetNumberOfStages(int, int, int);

    void                  InsertOpaqueImage(avtImage_p);

    void                  SetRayFunction(avtRayFunction *);
    void                  SetScreen(int, int);
    void                  SetSamplesPerRay(int);
    void                  SetBackgroundColor(const unsigned char [3]);
    void                  SetBackgroundMode(int mode);
    void                  SetGradientBackgroundColors(const double [3],
                                                      const double [3]);
    int                   GetSamplesPerRay(void)  { return samplesPerRay; };
    const int            *GetScreen(void)         { return screen; };

    void                  SetKernelBasedSampling(bool v)
                                    { kernelBasedSampling = v; };

  protected:
    avtViewInfo           view;

    int                   screen[2];
    int                   samplesPerRay;
    bool                  kernelBasedSampling;
    int                   backgroundMode;
    unsigned char         background[3];
    double                gradBG1[3];
    double                gradBG2[3];
    avtRayFunction       *rayfoo;
    
    avtImage_p            opaqueImage;

    virtual void          Execute(void);
    virtual avtPipelineSpecification_p
                          PerformRestriction(avtPipelineSpecification_p);
    static int            GetNumberOfDivisions(int, int, int);
};


#endif


