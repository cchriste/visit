// ************************************************************************* //
//                                avtDDF.h                                   //
// ************************************************************************* //

#ifndef AVT_DDF_H
#define AVT_DDF_H

#include <ddf_exports.h>

#include <vector>
#include <string>


class     vtkDataArray;
class     vtkDataSet;

class     avtBinningScheme;
class     avtDDFFunctionInfo;


// ****************************************************************************
//  Class: avtDDF
//
//  Purpose:
//      This class represents a derived data function, which allows for data
//      to be represented on new sampling spaces.  Further documentation about
//      DDFs can be found in the document: ...
//
//  Programmer: Hank Childs
//  Creation:   February 12, 2006
//
// ****************************************************************************

class AVTDDF_API avtDDF
{
  public:
                           avtDDF(avtDDFFunctionInfo *, float *);
    virtual               ~avtDDF();

    avtDDFFunctionInfo    *GetFunctionInfo(void) { return functionInfo; };
    vtkDataArray          *ApplyFunction(vtkDataSet *);
    vtkDataSet            *CreateGrid(void);

  protected: 
    avtDDFFunctionInfo    *functionInfo;
    float                 *vals;
};


#endif


