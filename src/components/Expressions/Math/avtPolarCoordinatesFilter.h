// ************************************************************************* //
//                         avtPolarCoordinatesFilter.h                       //
// ************************************************************************* //

#ifndef AVT_POLAR_COORDINATES_FILTER_H
#define AVT_POLAR_COORDINATES_FILTER_H


#include <avtSingleInputExpressionFilter.h>


class     vtkCell;
class     vtkCellDataToPointData;
class     vtkDataArray;
class     vtkDataSet;
class     vtkIdList;
class     vtkScalarData;


// ****************************************************************************
//  Class: avtPolarCoordinatesFilter
//
//  Purpose:
//      A filter that calculates the polar coordinates for each point.
//
//  Programmer: Hank Childs
//  Creation:   November 18, 2002
//
// ****************************************************************************

class EXPRESSION_API avtPolarCoordinatesFilter 
    : public avtSingleInputExpressionFilter
{
  public:
                              avtPolarCoordinatesFilter() {;};
    virtual                  ~avtPolarCoordinatesFilter() {;};

    virtual const char       *GetType(void)   
                                  { return "avtPolarCoordinatesFilter"; };
    virtual const char       *GetDescription(void)
                                  { return "Calculating polar coordinates."; };

  protected:
    virtual vtkDataArray     *DeriveVariable(vtkDataSet *);
    virtual bool              IsPointVariable(void)  { return true; };  
    virtual int               GetVariableDimension() { return 3; }
};


#endif


