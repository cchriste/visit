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

#ifndef PY_SPREADSHEETATTRIBUTES_H
#define PY_SPREADSHEETATTRIBUTES_H
#include <Python.h>
#include <SpreadsheetAttributes.h>

//
// Functions exposed to the VisIt module.
//
#define SPREADSHEETATTRIBUTES_NMETH 30
void           PySpreadsheetAttributes_StartUp(SpreadsheetAttributes *subj, void *data);
void           PySpreadsheetAttributes_CloseDown();
PyMethodDef *  PySpreadsheetAttributes_GetMethodTable(int *nMethods);
bool           PySpreadsheetAttributes_Check(PyObject *obj);
SpreadsheetAttributes *  PySpreadsheetAttributes_FromPyObject(PyObject *obj);
PyObject *     PySpreadsheetAttributes_New();
PyObject *     PySpreadsheetAttributes_Wrap(const SpreadsheetAttributes *attr);
void           PySpreadsheetAttributes_SetParent(PyObject *obj, PyObject *parent);
void           PySpreadsheetAttributes_SetDefaults(const SpreadsheetAttributes *atts);
std::string    PySpreadsheetAttributes_GetLogString();
std::string    PySpreadsheetAttributes_ToString(const SpreadsheetAttributes *, const char *);
PyObject *     PySpreadsheetAttributes_getattr(PyObject *self, char *name);
int            PySpreadsheetAttributes_setattr(PyObject *self, char *name, PyObject *args);
extern PyMethodDef PySpreadsheetAttributes_methods[SPREADSHEETATTRIBUTES_NMETH];

#endif

