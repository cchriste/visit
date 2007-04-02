#ifndef PY_CONTOURATTRIBUTES_H
#define PY_CONTOURATTRIBUTES_H
#include <Python.h>
#include <ContourAttributes.h>

//
// Functions exposed to the VisIt module.
//
void            PyContourAttributes_StartUp(ContourAttributes *subj, void *data);
void            PyContourAttributes_CloseDown();
PyMethodDef    *PyContourAttributes_GetMethodTable(int *nMethods);
bool            PyContourAttributes_Check(PyObject *obj);
ContourAttributes *PyContourAttributes_FromPyObject(PyObject *obj);
PyObject       *PyContourAttributes_NewPyObject();
PyObject       *PyContourAttributes_WrapPyObject(const ContourAttributes *attr);
void            PyContourAttributes_SetDefaults(const ContourAttributes *atts);
std::string     PyContourAttributes_GetLogString();
std::string     PyContourAttributes_ToString(const ContourAttributes *, const char *);

#endif

