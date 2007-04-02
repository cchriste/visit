#ifndef PY_MATERIALATTRIBUTES_H
#define PY_MATERIALATTRIBUTES_H
#include <Python.h>
#include <MaterialAttributes.h>

//
// Functions exposed to the VisIt module.
//
void            PyMaterialAttributes_StartUp(MaterialAttributes *subj, void *data);
void            PyMaterialAttributes_CloseDown();
PyMethodDef    *PyMaterialAttributes_GetMethodTable(int *nMethods);
bool            PyMaterialAttributes_Check(PyObject *obj);
MaterialAttributes *PyMaterialAttributes_FromPyObject(PyObject *obj);
PyObject       *PyMaterialAttributes_NewPyObject();
PyObject       *PyMaterialAttributes_WrapPyObject(const MaterialAttributes *attr);
void            PyMaterialAttributes_SetDefaults(const MaterialAttributes *atts);
std::string     PyMaterialAttributes_GetLogString();
std::string     PyMaterialAttributes_ToString(const MaterialAttributes *, const char *);

#endif

