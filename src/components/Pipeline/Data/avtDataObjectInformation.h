// ************************************************************************* //
//                           avtDataObjectInformation.h                      //
// ************************************************************************* //

#ifndef AVT_DATA_OBJECT_INFORMATION_H
#define AVT_DATA_OBJECT_INFORMATION_H
#include <pipeline_exports.h>
#include <ref_ptr.h>

#include <avtDataAttributes.h>
#include <avtDataValidity.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

class     avtDataObjectString;
class     avtDataObjectWriter;

// ****************************************************************************
//  Class: avtDataObjectInformation
//
//  Purpose:
//      An auxiliary class intended only to be used by avt data objects.  Its
//      purpose is to information about a data object in an encapsulated way
//      that translates across data objects.  It is divided into two classes,
//      one for describing the attributes of the dataset, the other for
//      describing the validity of the dataset.
//
//  Programmer: Hank Childs
//  Creation:   October 25, 2000
//  
//  Modifications:
//
//    Hank Childs, Sat Mar 24 15:14:42 PST 2001
//    Split class into two classes, blew away previous comments since they
//    now apply to avtDataAttributes and avtDataValidity.
//
// ****************************************************************************

class PIPELINE_API avtDataObjectInformation
{
  public:
                             avtDataObjectInformation();
    virtual                 ~avtDataObjectInformation();

    void                     Copy(const avtDataObjectInformation &);
 
    avtDataAttributes       &GetAttributes(void)   { return atts; };
    avtDataValidity         &GetValidity(void)     { return validity; };

    void                     Merge(const avtDataObjectInformation &);
    void                     ParallelMerge(const ref_ptr<avtDataObjectWriter>);

    void                     Write(avtDataObjectString &, 
                                   const avtDataObjectWriter *);
    int                      Read(char *);

  protected:
    avtDataAttributes        atts;
    avtDataValidity          validity;

  private:
    static int               objectCount;
    static void              InitializeMPIStuff();
    static void              FinalizeMPIStuff();
#ifdef PARALLEL
    static MPI_Op            mpiOpUnifyDobInfo;
#endif

};

#endif
