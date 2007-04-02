
// ************************************************************************* //
//                            avtAUXFileFileFormat.h                         //
// ************************************************************************* //

#ifndef AVT_AUX_FILE_FORMAT_H
#define AVT_AUX_FILE_FORMAT_H

#include <avtSTSDFileFormat.h>
#include <vtkFloatArray.h>

#include <string>

using std::string;


// ****************************************************************************
//  Class: avtAUXFileFileFormat
//
//  Purpose:
//      Reads in AUXFile files as a plugin to VisIt.
//
//  Programmer: miller -- generated by xml2avt
//  Creation:   Tue Mar 15 08:29:20 PDT 2005
//
// ****************************************************************************

class avtAUXFileFileFormat : public avtSTSDFileFormat
{
  public:
                       avtAUXFileFileFormat(const char *filename);
    virtual           ~avtAUXFileFileFormat() {;};

    virtual const char    *GetType(void)   { return "AUXFile"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(const char *);
    virtual vtkDataArray  *GetVar(const char *);

  protected:

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);

    string                 fileName;
    char *                 fileBuf;

    int                    sizeX, sizeY;

    vtkFloatArray         *fluence;

};


#endif
