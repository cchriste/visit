// ************************************************************************* //
//                            avtSimDBFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_SimDB_FILE_FORMAT_H
#define AVT_SimDB_FILE_FORMAT_H

#include <database_exports.h>

#include <avtSTMDFileFormat.h>

#include <vector>


// ****************************************************************************
//  Class: avtSimDBFileFormat
//
//  Purpose:
//      Reads in SimDB files as a plugin to VisIt.
//
//  Programmer: Jeremy Meredith
//  Creation:   March 25, 2004
//
// ****************************************************************************

class avtSimDBFileFormat : public avtSTMDFileFormat
{
  public:
                       avtSimDBFileFormat(const char *);
    virtual           ~avtSimDBFileFormat() {;};

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    // virtual void      &GetAuxiliaryData(const char *var, const char *type,
    //                                     int domain, void *args, 
    //                                     DestructorFunction &);
    //

    //
    // If you know the cycle number, overload this function.
    // Otherwise, VisIt will make up a reasonable one for you.
    //
    // virtual int         GetCyle(void);
    //

    virtual const char    *GetType(void)   { return "SimDB"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

  protected:
    string                 host;
    int                    port;

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);
};


#endif
