#ifndef PDBFILEOBJECT_H
#define PDBFILEOBJECT_H
#include <pdb.h>
#include <string>

//
// Enum for array type definitions.
//
typedef enum {NO_TYPE, CHAR_TYPE, INTEGER_TYPE, FLOAT_TYPE, DOUBLE_TYPE,
              LONG_TYPE, CHARARRAY_TYPE, INTEGERARRAY_TYPE, FLOATARRAY_TYPE,
              DOUBLEARRAY_TYPE, LONGARRAY_TYPE, OBJECT_TYPE} TypeEnum;

// ****************************************************************************
// Class: PDBFileObject
//
// Purpose:
//   This class encapsulates a PDB file and provides methods to read out
//   certain types of values.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Sep 16 09:00:58 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class PDBFileObject
{
public:
    PDBFileObject(const char *name);
    virtual ~PDBFileObject();

    bool IsOpen() const;
    bool Open();
    void Close();
    const std::string &GetName() const;
    PDBfile *filePointer();

    bool GetDouble(const char *name, double *val);
    bool GetDoubleArray(const char *name, double **val, int *nvals);
    bool GetInteger(const char *name, int *val);
    bool GetIntegerArray(const char *name, int **val, int *nvals);
    bool GetString(const char *name, char **str, int *len = 0);
    bool SymbolExists(const char *name);
    bool SymbolExists(const char *, TypeEnum *, int *, int **, int*);
    bool SymbolExists(const char *, TypeEnum *, std::string &, int *, int **, int*);
    void *ReadValues(const char *, TypeEnum *, int *, int **, int*, int=0);

private:
    bool AutoOpen();

    std::string  filename;
    PDBfile     *pdb;
};

//
// Functions to free memory.
//

template <class T>
void free_mem(T *);

void free_void_mem(void *, TypeEnum);

#endif
