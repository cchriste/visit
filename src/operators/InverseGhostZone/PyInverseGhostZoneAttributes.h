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

#ifndef PY_INVERSEGHOSTZONEATTRIBUTES_H
#define PY_INVERSEGHOSTZONEATTRIBUTES_H
#include <Python.h>
#include <InverseGhostZoneAttributes.h>

//
// Functions exposed to the VisIt module.
//
#define INVERSEGHOSTZONEATTRIBUTES_NMETH 16
void           PyInverseGhostZoneAttributes_StartUp(InverseGhostZoneAttributes *subj, void *data);
void           PyInverseGhostZoneAttributes_CloseDown();
PyMethodDef *  PyInverseGhostZoneAttributes_GetMethodTable(int *nMethods);
bool           PyInverseGhostZoneAttributes_Check(PyObject *obj);
InverseGhostZoneAttributes *  PyInverseGhostZoneAttributes_FromPyObject(PyObject *obj);
PyObject *     PyInverseGhostZoneAttributes_New();
PyObject *     PyInverseGhostZoneAttributes_Wrap(const InverseGhostZoneAttributes *attr);
void           PyInverseGhostZoneAttributes_SetParent(PyObject *obj, PyObject *parent);
void           PyInverseGhostZoneAttributes_SetDefaults(const InverseGhostZoneAttributes *atts);
std::string    PyInverseGhostZoneAttributes_GetLogString();
std::string    PyInverseGhostZoneAttributes_ToString(const InverseGhostZoneAttributes *, const char *);
PyObject *     PyInverseGhostZoneAttributes_getattr(PyObject *self, char *name);
int            PyInverseGhostZoneAttributes_setattr(PyObject *self, char *name, PyObject *args);
extern PyMethodDef PyInverseGhostZoneAttributes_methods[INVERSEGHOSTZONEATTRIBUTES_NMETH];

#endif

