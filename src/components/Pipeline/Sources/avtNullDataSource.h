// ************************************************************************* //
//                         avtNullDataSource.h                               //
// ************************************************************************* //

#ifndef AVT_NULL_DATA_SOURCE_H
#define AVT_NULL_DATA_SOURCE_H
#include <pipeline_exports.h>


#include <avtDataObjectSource.h>
#include <avtNullData.h>

// ****************************************************************************
//  Class: avtNullDataSource
//
//  Purpose:
//      A data object source whose data object is null data 
//
//  Programmer: Mark C. Miller 
//  Creation:   January 8, 2003 
//
// ****************************************************************************

class PIPELINE_API avtNullDataSource : virtual public avtDataObjectSource
{
  public:
                                avtNullDataSource();
    virtual                    ~avtNullDataSource() {;};

    avtDataObject_p             GetOutput(void);
    avtNullData_p               GetTypedOutput(void) { return nullData; };

  protected:
    avtNullData_p               nullData;

    void                        SetOutput(avtNullData_p);

};


#endif


