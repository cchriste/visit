#ifndef PP_Z_STSD_FILE_FORMAT_H
#define PP_Z_STSD_FILE_FORMAT_H
#include <avtSTSDFileFormat.h>
#include <PP_ZFileReader.h>

class avtFileFormatInterface;

// ****************************************************************************
// Class: PP_Z_STSD_FileFormat
//
// Purpose:
//   This file format uses
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Sep 16 15:27:25 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class PP_Z_STSD_FileFormat : public avtSTSDFileFormat
{
public:
    static avtFileFormatInterface *CreateInterface(PDBFileObject *pdb,
         const char *const *filenames, int nList, int nBlocks);

    PP_Z_STSD_FileFormat(const char *filename);
    PP_Z_STSD_FileFormat(PDBFileObject *p);
    virtual ~PP_Z_STSD_FileFormat();

    // Mimic PDBReader interface.
    bool Identify();
    void SetOwnsPDBFile(bool);

    // Methods overrides for an STSD file format.
    virtual const char   *GetType();
    virtual void          PopulateDatabaseMetaData(avtDatabaseMetaData *);
    virtual int           GetCycle(void);
    virtual vtkDataSet   *GetMesh(const char *var);
    virtual vtkDataArray *GetVar(const char *var);
    virtual void         *GetAuxiliaryData(const char *var, const char *type,
                                           void *args, DestructorFunction &);
protected:
    PP_ZFileReader        reader;
};

#endif
