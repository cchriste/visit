// ************************************************************************* //
//                            avtNullDataWriter.h                            //
// ************************************************************************* //

#ifndef AVT_NULL_DATA_WRITER_H
#define AVT_NULL_DATA_WRITER_H

#include <pipeline_exports.h>

#include <avtOriginatingNullDataSink.h>
#include <avtDataObjectWriter.h>


class     avtDataObjectString;


// ****************************************************************************
//  Class: avtNullDataWriter
//
//  Purpose:
//      A class which takes as input an avtNullData and can serialize it.
//
//  Programmer: Mark C. Miller 
//  Creation:   January 7, 2003 
//
//  Modifications:
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined constructor and destructor definitions to .C files
//    because certain compilers have problems with them.
//
// ****************************************************************************

class PIPELINE_API avtNullDataWriter : public avtOriginatingNullDataSink,
                       public avtDataObjectWriter
{
  public:
                       avtNullDataWriter();
    virtual           ~avtNullDataWriter();

    //virtual bool       MustMergeParallelStreams(void) { return true; };

  protected:
    void               DataObjectWrite(avtDataObjectString &);
};


#endif


