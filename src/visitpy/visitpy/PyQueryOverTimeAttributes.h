/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef PY_QUERYOVERTIMEATTRIBUTES_H
#define PY_QUERYOVERTIMEATTRIBUTES_H
#include <Python.h>
#include <QueryOverTimeAttributes.h>
#include <visitpy_exports.h>

//
// Functions exposed to the VisIt module.
//
void VISITPY_API           PyQueryOverTimeAttributes_StartUp(QueryOverTimeAttributes *subj, void *data);
void VISITPY_API           PyQueryOverTimeAttributes_CloseDown();
VISITPY_API PyMethodDef *   PyQueryOverTimeAttributes_GetMethodTable(int *nMethods);
bool VISITPY_API           PyQueryOverTimeAttributes_Check(PyObject *obj);
VISITPY_API QueryOverTimeAttributes *  PyQueryOverTimeAttributes_FromPyObject(PyObject *obj);
VISITPY_API PyObject *      PyQueryOverTimeAttributes_New();
VISITPY_API PyObject *      PyQueryOverTimeAttributes_Wrap(const QueryOverTimeAttributes *attr);
void VISITPY_API           PyQueryOverTimeAttributes_SetParent(PyObject *obj, PyObject *parent);
void VISITPY_API           PyQueryOverTimeAttributes_SetDefaults(const QueryOverTimeAttributes *atts);
std::string VISITPY_API    PyQueryOverTimeAttributes_GetLogString();
std::string VISITPY_API    PyQueryOverTimeAttributes_ToString(const QueryOverTimeAttributes *, const char *);

#endif

