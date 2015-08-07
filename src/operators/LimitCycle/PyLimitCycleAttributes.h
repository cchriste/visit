/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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

#ifndef PY_LIMITCYCLEATTRIBUTES_H
#define PY_LIMITCYCLEATTRIBUTES_H
#include <Python.h>
#include <LimitCycleAttributes.h>

//
// Functions exposed to the VisIt module.
//
#define LIMITCYCLEATTRIBUTES_NMETH 120
void           PyLimitCycleAttributes_StartUp(LimitCycleAttributes *subj, void *data);
void           PyLimitCycleAttributes_CloseDown();
PyMethodDef *  PyLimitCycleAttributes_GetMethodTable(int *nMethods);
bool           PyLimitCycleAttributes_Check(PyObject *obj);
LimitCycleAttributes *  PyLimitCycleAttributes_FromPyObject(PyObject *obj);
PyObject *     PyLimitCycleAttributes_New();
PyObject *     PyLimitCycleAttributes_Wrap(const LimitCycleAttributes *attr);
void           PyLimitCycleAttributes_SetParent(PyObject *obj, PyObject *parent);
void           PyLimitCycleAttributes_SetDefaults(const LimitCycleAttributes *atts);
std::string    PyLimitCycleAttributes_GetLogString();
std::string    PyLimitCycleAttributes_ToString(const LimitCycleAttributes *, const char *);
PyObject *     PyLimitCycleAttributes_getattr(PyObject *self, char *name);
int            PyLimitCycleAttributes_setattr(PyObject *self, char *name, PyObject *args);
extern PyMethodDef PyLimitCycleAttributes_methods[LIMITCYCLEATTRIBUTES_NMETH];

#endif

