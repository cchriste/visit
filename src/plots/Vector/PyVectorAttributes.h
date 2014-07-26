/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef PY_VECTORATTRIBUTES_H
#define PY_VECTORATTRIBUTES_H
#include <Python.h>
#include <VectorAttributes.h>

//
// Functions exposed to the VisIt module.
//
#define VECTORATTRIBUTES_NMETH 56
void           PyVectorAttributes_StartUp(VectorAttributes *subj, void *data);
void           PyVectorAttributes_CloseDown();
PyMethodDef *  PyVectorAttributes_GetMethodTable(int *nMethods);
bool           PyVectorAttributes_Check(PyObject *obj);
VectorAttributes *  PyVectorAttributes_FromPyObject(PyObject *obj);
PyObject *     PyVectorAttributes_New();
PyObject *     PyVectorAttributes_Wrap(const VectorAttributes *attr);
void           PyVectorAttributes_SetParent(PyObject *obj, PyObject *parent);
void           PyVectorAttributes_SetDefaults(const VectorAttributes *atts);
std::string    PyVectorAttributes_GetLogString();
std::string    PyVectorAttributes_ToString(const VectorAttributes *, const char *);
PyObject *     PyVectorAttributes_getattr(PyObject *self, char *name);
int            PyVectorAttributes_setattr(PyObject *self, char *name, PyObject *args);
extern PyMethodDef PyVectorAttributes_methods[VECTORATTRIBUTES_NMETH];

#endif

