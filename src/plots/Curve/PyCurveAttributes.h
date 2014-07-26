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

#ifndef PY_CURVEATTRIBUTES_H
#define PY_CURVEATTRIBUTES_H
#include <Python.h>
#include <CurveAttributes.h>

//
// Functions exposed to the VisIt module.
//
#define CURVEATTRIBUTES_NMETH 58
void           PyCurveAttributes_StartUp(CurveAttributes *subj, void *data);
void           PyCurveAttributes_CloseDown();
PyMethodDef *  PyCurveAttributes_GetMethodTable(int *nMethods);
bool           PyCurveAttributes_Check(PyObject *obj);
CurveAttributes *  PyCurveAttributes_FromPyObject(PyObject *obj);
PyObject *     PyCurveAttributes_New();
PyObject *     PyCurveAttributes_Wrap(const CurveAttributes *attr);
void           PyCurveAttributes_SetParent(PyObject *obj, PyObject *parent);
void           PyCurveAttributes_SetDefaults(const CurveAttributes *atts);
std::string    PyCurveAttributes_GetLogString();
std::string    PyCurveAttributes_ToString(const CurveAttributes *, const char *);
PyObject *     PyCurveAttributes_getattr(PyObject *self, char *name);
int            PyCurveAttributes_setattr(PyObject *self, char *name, PyObject *args);
extern PyMethodDef PyCurveAttributes_methods[CURVEATTRIBUTES_NMETH];

#endif

