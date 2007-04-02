// ************************************************************************* //
//                             avtSTSDFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_STSD_FILE_FORMAT_H
#define AVT_STSD_FILE_FORMAT_H

#include <database_exports.h>

#include <avtFileFormat.h>

#include <void_ref_ptr.h>

class     vtkDataSet;
class     vtkDataArray;

class     avtIOInformation;


// ****************************************************************************
//  Class: avtSTSDFileFormat
//
//  Purpose:
//      This is an implementation of avtFileFormat for file formats that have
//      a Single Timestep and a Single Domain.
//
//  Programmer: Hank Childs
//  Creation:   February 22, 2001
//
//  Modifications:
//
//    Hank Childs, Mon Mar 11 08:48:41 PST 2002
//    Removed SetDatabaseMetaData since it is defined in the base class.
//
//    Kathleen Bonnell, Mon Mar 18 17:22:30 PST 2002  
//    vtkScalars and vtkVectors have been deprecated in VTK 4.0, 
//    use vtkDataArray instead.
//
//    Kathleen Bonnell, Mon Mar 18 17:22:30 PST 2002  
//    vtkScalars and vtkVectors have been deprecated in VTK 4.0, 
//    use vtkDataArray instead.
//
//    Brad Whitlock, Mon Oct 13 14:42:42 PST 2003
//    Added a method to return the time. The default implementation just
//    returns zero but it allows derived types to perform other actions
//    to return the right time. I also added method to indicate whether the
//    returned cycle/time can be believed.
//
//    Mark C. Miller, Mon Feb 23 20:38:47 PST 2004
//    Added method, ActivateTimestep
//
//    Mark C. Miller, Tue Mar 16 14:28:42 PST 2004
//    Added method, PopulateIOInformation
//
//    Mark C. Miller, Tue May 17 18:48:38 PDT 2005
//    Removed ReturnsValidCycle/Time GetCycle/Time. Added SetDatabaseMetaData
//
// ****************************************************************************

class DATABASE_API avtSTSDFileFormat : public avtFileFormat
{
  public:
                           avtSTSDFileFormat(const char *);
    virtual               ~avtSTSDFileFormat();

    void                   SetTimestep(int ts) { timestep = ts; };
    void                   SetDomain(int dom) { domain = dom; };
    virtual void          *GetAuxiliaryData(const char *var, const char *type,
                                            void *args, DestructorFunction &);

    virtual const char    *GetFilename(void) { return filename; };

    virtual vtkDataSet    *GetMesh(const char *) = 0;
    virtual vtkDataArray  *GetVar(const char *) = 0;
    virtual vtkDataArray  *GetVectorVar(const char *);

    virtual void           ActivateTimestep(void)
                               { avtFileFormat::ActivateTimestep(); };
    virtual void           PopulateIOInformation(avtIOInformation& ioInfo)
                               { avtFileFormat::PopulateIOInformation(ioInfo); };
    virtual void           SetDatabaseMetaData(avtDatabaseMetaData *md)
                               { metadata = md; PopulateDatabaseMetaData(metadata); };

  protected:
    char                  *filename;
    int                    domain;
    int                    timestep;

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData*) = 0;

};


#endif


