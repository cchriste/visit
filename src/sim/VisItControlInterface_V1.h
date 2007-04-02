#ifndef VISIT_CONTROL_INTERFACE_V1_H
#define VISIT_CONTROL_INTERFACE_V1_H

/* ****************************************************************************
//  File:  VisItControlInterfave_V1.h
//
//  Purpose:
//    Abstraction of VisIt Engine wrapper library.  Handles the
//    grunt work of actually connecting to visit that must be done
//    outside of the VisItEngine DLL.
//
//  Programmer:  Jeremy Meredith
//  Creation:    April  4, 2005
//
// ***************************************************************************/
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
                                          const char *inputfile);
int   VisItDetectInput(int blocking, int consoledesc);
int   VisItAttemptToCompleteConnection(void);
void  VisItSetSlaveProcessCallback(void(*)());
void  VisItSetCommandCallback(void(*)(const char*,int,float,const char*));
int   VisItProcessEngineCommand(void);
void  VisItTimeStepChanged(void);
void  VisItDisconnect(void);
char *VisItGetLastError();

#ifdef __cplusplus
}
#endif

#endif
