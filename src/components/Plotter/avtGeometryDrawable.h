// ************************************************************************* //
//                           avtGeometryDrawable.h                           //
// ************************************************************************* //

#ifndef AVT_GEOMETRY_DRAWABLE_H
#define AVT_GEOMETRY_DRAWABLE_H

#include <plotter_exports.h>

#include <avtDrawable.h>


class     vtkActor;
class     vtkDataObjectCollection;

class     avtMapper;

class     ColorAttribute;

// ****************************************************************************
//  Class: avtGeometryDrawable
//
//  Purpose:
//      A concrete type of avtDrawable, this allows for adding and removing
//      actors that come from geometry to/from a renderer.
//
//  Programmer: Hank Childs
//  Creation:   December 21, 2000
//
//  Modifications:
//
//    Hank Childs, Mon Mar 12 16:54:39 PST 2001
//    Allowed for geometry to be shifted by a vector.
//
//    Kathleen Bonnell, Thu Mar 15 10:51:59 PST 2001
//    Added member actorsVisibility to store visible state of individual actors.
//
//    Hank Childs, Sun Jul  7 14:16:38 PDT 2002
//    Add support for transparency.
//
//    Kathleen Bonnell, Fri Jul 12 16:21:37 PDT 2002
//    Added method ScaleByVector.
//
//    Kathleen Bonnell, Fri Jul 19 08:39:04 PDT 2002 
//    Added method UpdateScaleFactor.
//
//    Kathleen Bonnell, Tue Aug 13 15:15:37 PDT 2002   
//    Added methods in support of lighting. 
//
//    Brad Whitlock, Mon Sep 23 15:57:46 PST 2002
//    Added a method to set the surface representation. I added another method
//    for setting the drawable's immediate mode rendering flag.
//
//    Jeremy Meredith, Fri Nov 14 11:10:53 PST 2003
//    Added a method to set the specular properties.
//
//    Mark C. Miller, Tue May 11 20:21:24 PDT 2004
//    Removed SetExternallyRenderedImagesActor method
//
//    Kathleen Bonnell, Mon Sep 27 10:21:15 PDT 2004 
//    Added MakePickable and MakeUnPickable. 
//
//    Kathleen Bonnell, Mon Jun 27 15:19:14 PDT 2005 
//    Added GetZPosition.
//
// ****************************************************************************

class PLOTTER_API avtGeometryDrawable : public avtDrawable
{
  public:
                                avtGeometryDrawable(int, vtkActor **);
    virtual                    ~avtGeometryDrawable();

    void                        SetMapper(avtMapper *);

    virtual bool                Interactive(void)  { return true; };

    virtual void                Add(vtkRenderer *);
    virtual void                Remove(vtkRenderer *);

    virtual void                VisibilityOn(void);
    virtual void                VisibilityOff(void);
    virtual int                 SetTransparencyActor(avtTransparencyActor *);

    virtual void                ShiftByVector(const double [3]);
    virtual void                ScaleByVector(const double [3]);
    virtual void                UpdateScaleFactor(void);

    virtual void                TurnLightingOn(void);
    virtual void                TurnLightingOff(void);
    virtual void                SetAmbientCoefficient(const double);

    virtual void                SetSurfaceRepresentation(int rep);
    virtual void                SetImmediateModeRendering(bool val);

    virtual void                SetSpecularProperties(bool,double,double,
                                                      const ColorAttribute&);

    virtual avtDataObject_p     GetDataObject(void);

    virtual void                MakePickable(void);
    virtual void                MakeUnPickable(void);
    virtual double              GetZPosition(void);

  protected:
    int                         nActors;
    vtkActor                  **actors;
    vtkRenderer                *renderer;
    int                        *actorsVisibility;
    avtMapper                  *mapper;
};


#endif


