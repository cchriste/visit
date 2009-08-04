/*****************************************************************************
*
* Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400124
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

#ifndef VISIT_CONTROL_INTERFACE_V2_H
#define VISIT_CONTROL_INTERFACE_V2_H
#include <VisItInterfaceTypes_V2.h>

/*****************************************************************************
 *  File:  VisItControlInterface_V2.h
 *
 *  Purpose:
 *    Abstraction of VisIt Engine wrapper library.  Handles the
 *    grunt work of actually connecting to visit that must be done
 *    outside of the VisItEngine DLL.
 *
 *  Programmer:  Jeremy Meredith
 *  Creation:    April  4, 2005
 *
 *  Modifications:
 *    Shelly Prevost, Wed Jan 25 08:52:18 PST 2006
 *    Added a guifile argument to VisItInitializeSocketAndDumpSimFile.
 *
 *    Brad Whitlock, Thu Jan 25 14:53:11 PST 2007
 *    Added VisItUpdatePlots and VisItExecuteCommand.
 *
 *    Brad Whitlock, Mon Feb  9 10:04:49 PST 2009
 *    Added new functions, upped version to "V2".
 *
 *****************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

void  VisItSetBroadcastIntFunction(int (*)(int *, int));
void  VisItSetBroadcastStringFunction(int (*)(char *, int, int));
void  VisItSetParallel(int);
void  VisItSetParallelRank(int);

void  VisItSetDirectory(char *);
void  VisItSetOptions(char *);
int   VisItSetupEnvironment(void);
int   VisItInitializeSocketAndDumpSimFile(const char *name,
                                          const char *comment,
                                          const char *path,
                                          const char *inputfile,
                                          const char *guifile,
                                          const char *absoluteFilename);
int   VisItDetectInput(int blocking, int consoledesc);
int   VisItAttemptToCompleteConnection(void);
void  VisItSetSlaveProcessCallback(void(*)(void));
void  VisItSetCommandCallback(void(*)(const char*,const char*,void*), void*);
int   VisItProcessEngineCommand(void);
void  VisItTimeStepChanged(void);
void  VisItUpdatePlots(void);
void  VisItExecuteCommand(const char *);
void  VisItDisconnect(void);
int   VisItIsConnected(void);
char *VisItGetLastError(void);
int   VisItSynchronize(void);
void  VisItEnableSynchronize(int);

void  VisItDebug1(const char *format, ...);
void  VisItDebug2(const char *format, ...);
void  VisItDebug3(const char *format, ...);
void  VisItDebug4(const char *format, ...);
void  VisItDebug5(const char *format, ...);

void  VisItOpenTraceFile(const char *);
void  VisItCloseTraceFile(void);

int VisItSaveWindow(const char *filename, int width, int height, int format);

/* This would come from the data interface */
/*typedef struct VisIt_MeshData_t* VisIt_MeshData;*/
#include <VisItDataInterface_V2.h>

/* Functions that install data access callback functions */

/* This is needed to do collective communication before we call the other callbacks. */
int VisItSetActivateTimestep(int (*cb)(void *), void *cbdata);
int VisItSetGetMetaData(int (*cb)(VisIt_SimulationMetaData *, void *), void *cbdata);
int VisItSetGetMesh(int (*cb)(int, const char *, VisIt_MeshData *, void *), void *cbdata);
int VisItSetGetMaterial(int (*cb)(int, const char *, VisIt_MaterialData *, void *), void *cbdata);
int VisItSetGetSpecies(int (*cb)(int, const char *, VisIt_SpeciesData *, void *), void *cbdata);
int VisItSetGetVariable(int (*cb)(int, const char *, VisIt_VariableData *, void *), void *cbdata);
int VisItSetGetMixedVariable(int (*cb)(int, const char *, VisIt_MixedVariableData *, void *), void *cbdata);
int VisItSetGetCurve(int (*cb)(const char *, VisIt_CurveData *, void *), void *cbdata);
int VisItSetGetDomainList(int (*cb)(VisIt_DomainList *, void *), void *cbdata);

/* This is needed for VisIt to create ghost zones in between the domains. */
int VisItSetGetDomainBoundaries(int (*cb)(const char *, visit_handle, void *), void *cbdata);

/* This is needed to tell VisIt how AMR patches are nested. */
int VisItSetGetDomainNesting(int (*cb)(const char *, visit_handle, void *), void *cbdata);

/* Functions that install data writer callback functions */
int VisItSetWriteBegin(int (*cb)(void *, const char *), void *cbdata);
int VisItSetWriteEnd(int (*cb)(void *, const char *), void *cbdata);
int VisItSetWriteMesh(int (*cb)(void *, const char *, int, const VisIt_MeshData *, const VisIt_MeshMetaData *), void *cbdata);
int VisItSetWriteVariable(int (*cb)(void *, const char *, const char *, int, int, void *, int, int, const VisIt_VariableMetaData *), void *cbdata);

#ifdef __cplusplus
}
#endif

#endif
