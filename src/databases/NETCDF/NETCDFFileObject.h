#ifndef NETCDF_FILE_OBJECT_H
#define NETCDF_FILE_OBJECT_H
#include <string>
#include <visitstream.h>

typedef enum {NO_TYPE, CHARARRAY_TYPE, UCHARARRAY_TYPE, SHORTARRAY_TYPE,
              INTEGERARRAY_TYPE, FLOATARRAY_TYPE,
              DOUBLEARRAY_TYPE, LONGARRAY_TYPE} TypeEnum;

// ****************************************************************************
// Class: NETCDFFileObject
//
// Purpose:
//   Abstract the NETCDF file a little.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri Aug 12 10:07:08 PDT 2005
//
// Modifications:
//   
// ****************************************************************************

class NETCDFFileObject
{
public:
    NETCDFFileObject(const char *name); 
    virtual ~NETCDFFileObject();

    bool IsOpen() const;
    bool Open();
    void Close();
    const std::string &GetName() const;
    int  GetFileHandle();

    bool ReadAttribute(const char *attname, TypeEnum *type, int *ndims,
                       int **dims, void **value);
    bool ReadAttribute(const char *varname, const char *attname,
                       TypeEnum *type, int *ndims, int **dims, void **value);
    // Convenience functions
    bool ReadStringAttribute(const char *varname, const char *attname,
                             std::string &attval);
    bool ReadStringAttribute(const char *attname, std::string &attval);

    bool InqVariable(const char *varname, TypeEnum *, int *ndims, int **dims);
    bool ReadVariable(const char *varname, TypeEnum *type, int *ndims,
                      int **dims, void **values);
    bool ReadVariableInto(const char *varname, TypeEnum t, void *arr);
    bool ReadVariableIntoAsFloat(const char *varname, float *arr);

    bool GetVarId(const char *name, int *varid);

    void HandleError(int status) const;
    void PrintFileContents(ostream &os);
private:
    bool AutoOpen();

    std::string filename;
    int         ncid;
};

//
// Functions to free memory.
//

template <class T>
void free_mem(T *);

void free_void_mem(void *, TypeEnum);

#endif
