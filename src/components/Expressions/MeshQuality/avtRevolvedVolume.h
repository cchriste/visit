// ************************************************************************* //
//                          avtRevolvedVolume.h                              //
// ************************************************************************* //

#ifndef AVT_REVOLVED_VOLUME_H
#define AVT_REVOLVED_VOLUME_H


#include <avtSingleInputExpressionFilter.h>

class     vtkCell;


// ****************************************************************************
//  Class: avtRevolvedVolume
//
//  Purpose:
//      Calculates the volume a 2D polygon would occupy if it were revolved
//      around an axis.
//
//  Programmer: Hank Childs
//  Creation:   September 8, 2002
//
// ****************************************************************************

class EXPRESSION_API avtRevolvedVolume : public avtSingleInputExpressionFilter
{
  public:
                                avtRevolvedVolume();

    virtual const char         *GetType(void) { return "avtRevolvedVolume"; };
    virtual const char         *GetDescription(void)
                                    { return "Calculating revolved volume"; };
    
  protected:
    bool                        haveIssuedWarning;

    virtual vtkDataArray       *DeriveVariable(vtkDataSet *);
    virtual void                PreExecute(void);

    virtual bool                IsPointVariable(void)  { return false; };

    double                      GetZoneVolume(vtkCell *);
    double                      GetTriangleVolume(double [3], double [3]);
    double                      RevolveLineSegment(double [2], double [2],
                                                   double *);
};


#endif


