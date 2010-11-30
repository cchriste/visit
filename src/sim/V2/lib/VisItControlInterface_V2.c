/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
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

#include "VisItControlInterface_V2.h"
#include "SimV2Tracing.h"
#include "DeclareDataCallbacks.h"

#ifdef _WIN32
#include <winsock2.h>
#include <direct.h>
#else
#include <dlfcn.h>
#include <netdb.h>
#include <unistd.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <pwd.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#endif
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>

#include <visit-config.h> /* For HAVE_SOCKLEN_T */

#ifdef _WIN32
#define SNPRINTF _snprintf
#else
#define SNPRINTF sprintf
#endif

/* ****************************************************************************
 *  File:  VisItControlInterface.c
 *
 *  Purpose:
 *    Abstraction of VisIt Engine wrapper library.  Handles the
 *    grunt work of actually connecting to visit that must be done
 *    outside of the VisItEngine DLL, such as:
 *       1) setting up a listen socket
 *       2) writing a .sim file
 *       3) opening the VisItEngine .so and retrieving the functions from it
 *       4) accepting an incoming socket connection
 *       5) removing the .sim file when the program exits
 *
 *  Programmer:  Jeremy Meredith
 *  Creation:    May 5, 2005
 *
 *****************************************************************************/

#define INITIAL_PORT_NUMBER 5609

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

static int BroadcastInt(int *value, int sender);


/* VisIt Engine Library function pointers */
typedef struct
{
    void *(*get_engine)(void);
    int   (*get_descriptor)(void*);
    int   (*process_input)(void*);
    int   (*initialize)(void*,int,char**);
    int   (*connect_viewer)(void*,int,char**);
    void  (*time_step_changed)(void*);
    void  (*execute_command)(void*,const char*);
    void  (*disconnect)();
    void  (*set_slave_process_callback)(void(*)());
    void  (*set_command_callback)(void*,void(*)(const char*,const char*,void*),void*);
    int   (*save_window)(void*,const char *, int, int, int);
    void  (*debug_logs)(int,const char *);
} control_callback_t;

#define STRUCT_MEMBER(F,T)  void (*set_##F)(T, void*);

typedef struct
{
    DECLARE_DATA_CALLBACKS(STRUCT_MEMBER, NO_ARGS)
} data_callback_t;

typedef struct
{
    control_callback_t control;
    data_callback_t    data;
} visit_callback_t;

typedef struct
{
    int    id;
    void (*cb)(void*);
    void  *cbdata;
} SyncCallback;

static visit_callback_t *callbacks = NULL;

/* Internal Variables */
#ifdef _WIN32
#define VISIT_SOCKET         SOCKET
#define VISIT_INVALID_SOCKET INVALID_SOCKET

static HMODULE     dl_handle;
static SOCKET      listenSocket = VISIT_INVALID_SOCKET;
static SOCKET      engineSocket = VISIT_INVALID_SOCKET;
#else
#define VISIT_SOCKET         int
#define VISIT_INVALID_SOCKET -1

static void       *dl_handle = NULL;
static int         listenSocket = VISIT_INVALID_SOCKET;
static int         engineSocket = VISIT_INVALID_SOCKET;
#endif

static int       (*BroadcastInt_internal)(int *value, int sender) = NULL;
static int       (*BroadcastString_internal)(char *str, int len, int sender) = NULL;
static char       *visit_directory = NULL;
static char       *visit_options = NULL;
static void       *engine = NULL;
static int         engine_argc = 0;
static char      **engine_argv;
static char        simulationFileName[1024];
static char        securityKey[17];
static char        localhost[256];
static int         listenPort = -1;
struct sockaddr_in listenSockAddr;
static int         isParallel = FALSE;
static int         parallelRank = 0;
static char        lastError[1024] = "";
static int           visit_sync_enabled = 1;
static SyncCallback *visit_sync_callbacks = NULL;
static int           visit_sync_callbacks_size = 0;
static int           visit_sync_id = 1;
static void        (*visit_command_callback)(const char*,const char*,void*) = NULL;
static void         *visit_command_callback_data = NULL;


/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
                               Internal Functions
 *******************************************************************************
 *******************************************************************************
 ******************************************************************************/
#ifdef _WIN32
/*******************************************************************************
* UTILITY FUNCTIONS
*******************************************************************************/
static int
ReadKeyFromRoot(HKEY which_root, const char *ver, const char *key,
    char **keyval)
{
    int  readSuccess = 0;
    char regkey[100];
    HKEY hkey;

    /* Try and read the key from the system registry. */
    sprintf(regkey, "VISIT%s", ver);
    *keyval = (char *)malloc(500);
    if(RegOpenKeyEx(which_root, regkey, 0, KEY_QUERY_VALUE, &hkey) == ERROR_SUCCESS)
    {
        DWORD keyType, strSize = 500;
        if(RegQueryValueEx(hkey, key, NULL, &keyType,
           (unsigned char *)*keyval, &strSize) == ERROR_SUCCESS)
        {
            readSuccess = 1;
        }

        RegCloseKey(hkey);
    }

    return readSuccess;
}

static int
ReadKey(const char *ver, const char *key, char **keyval)
{
    int retval = 0;

    if((retval = ReadKeyFromRoot(HKEY_CLASSES_ROOT, ver, key, keyval)) == 0)
        retval = ReadKeyFromRoot(HKEY_CURRENT_USER, ver, key, keyval);
    
    return retval;     
}

static void
GetVisItDirectory(char *visitdir, int maxlen)
{
    char *visitpath = 0;
    char visitversion[10];
    int major, minor, patch;
    int haveVISITHOME = 0;

    /* Iterate through a few probable versions and keep the most up to date one */
    for(major = 2; major <= 4; ++major)
        for(minor = (major==2?2:0); minor <= 12; ++minor)
            for(patch = 0; patch < 5; ++patch)
            {
                char curversion[10], *path = NULL;
                SNPRINTF(curversion, 10, "%d.%d.%d", major, minor, patch);
                if(ReadKey(curversion, "VISITHOME", &path))
                {
                    strcpy(visitversion, curversion);
                    if(visitpath != NULL)
                        free(visitpath);
                    visitpath = path;
                    haveVISITHOME = 1;
                }
            }
    if(haveVISITHOME)
    {
        SNPRINTF(visitdir, maxlen, "%s", visitpath);
        free(visitpath);
    }
    else if(visit_directory != NULL)
    {
        SNPRINTF(visitdir, maxlen, "%s", visit_directory);
    }
    else
    {
        visitdir[0] = '\0';
    }
}
#else
/* UNIX */
/*******************************************************************************
*
* Name: ReadEnvironmentFromCommand
*
* Purpose: Read the output of "visit -env" using the specified VisIt command.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Eric Brugger, Thu Sep 14 12:53:46 PDT 2006
*   Changed the routine to read at most ENV_BUF_SIZE bytes of output from
*   the execution of the command.
*
*   Brad Whitlock, Fri Jul 25 12:11:16 PDT 2008
*   Changed some types to remove warnings. Added trace information.
*
*******************************************************************************/
#define ENV_BUF_SIZE 10000

static int ReadEnvironmentFromCommand(const char *visitpath, char *output)
{
   /* VisIt will tell us what variables to set. */
   /* (redirect stderr so it won't complain if it can't find visit) */
   ssize_t n;
   size_t  lbuf;
   char command[1024];
   char *ptr;
   FILE *file;

   LIBSIM_API_ENTER1(ReadEnvironmentFromCommand, "visitpath=%s", visitpath);

#ifdef VISIT_COMPILER
#define STR(s) STR2(s)
#define STR2(s) #s
   SNPRINTF(command, 1024, "%s -compiler %s %s -env -engine 2>/dev/null",
           visitpath, STR(VISIT_COMPILER), visit_options ? visit_options : "");
#else
   SNPRINTF(command, 1024, "%s %s -env -engine 2>/dev/null",
           visitpath, visit_options ? visit_options : "");
#endif

   LIBSIM_MESSAGE1("command=%s", command);

   file = popen(command, "r");
   ptr = output;
   lbuf = ENV_BUF_SIZE;
   while ((n = read(fileno(file), (void*)ptr, lbuf)) > 0)
   {
      ptr += n;
      lbuf -= n;
   }
   *ptr = '\0';

   LIBSIM_MESSAGE1("Output=%s", output);

   LIBSIM_API_LEAVE1(ReadEnvironmentFromCommand, "return %d", (int)(ptr-output));
   return (ptr - output);
}

#endif

static void
VisItMkdir(const char *dir, int permissions)
{
#ifdef _WIN32
    _mkdir(dir);
#else
    mkdir(dir, permissions);
#endif
}



/******************************************************************************
*
* Name: visit_get_sync
*
* Purpose: This function returns a slot in the sync table, creating one if needed.
*
* Programmer: Brad Whitlock
* Date:       Thu Mar 26 13:28:11 PDT 2009
*
* Modifications:
*  Cyrus Harrison, Wed Nov 18 11:56:37 PST 2009
*  Declare "*ptr" at top of function for strict c compile.
*
******************************************************************************/

static SyncCallback *
visit_get_sync(void)
{
    SyncCallback *retval = NULL;
    SyncCallback *ptr = NULL;

    if(visit_sync_callbacks == NULL)
    {
        visit_sync_callbacks = (SyncCallback *)calloc(20, sizeof(SyncCallback));
        visit_sync_callbacks_size = 20;
        retval = visit_sync_callbacks;
    }
    else
    {
        /* return first empty slot. */
        int i;
        for(i = 0; i < visit_sync_callbacks_size; ++i)
        {
            if(visit_sync_callbacks[i].id == 0)
                return &visit_sync_callbacks[i];
        }

        /* resize */
        ptr = (SyncCallback *)calloc(visit_sync_callbacks_size+20, sizeof(SyncCallback));
        memcpy((void*)ptr, (void*)visit_sync_callbacks, visit_sync_callbacks_size * sizeof(SyncCallback));
        free(visit_sync_callbacks);
        visit_sync_callbacks = ptr;
        retval = &visit_sync_callbacks[visit_sync_callbacks_size];
        visit_sync_callbacks_size += 20;
    }
    return retval;
}

/******************************************************************************
*
* Name: visit_get_sync2
*
* Purpose: This function returns the slot with a given id.
*
* Programmer: Brad Whitlock
* Date:       Thu Mar 26 13:28:11 PDT 2009
*
* Modifications:
*
******************************************************************************/

static SyncCallback *
visit_get_sync2(int id)
{
    int i;
    for(i = 0; i < visit_sync_callbacks_size; ++i)
    {
        if(visit_sync_callbacks[i].id == id)
            return &visit_sync_callbacks[i];
    }
    return NULL;
}

/******************************************************************************
*
* Name: visit_add_sync
*
* Purpose: This function adds a callback function to the sync table and sends
*          a sync message to VisIt. When we get that message back, via the 
*          command callback mechanism, we execute the callback.
*
* Programmer: Brad Whitlock
* Date:       Thu Mar 26 13:28:11 PDT 2009
*
* Modifications:
*  Cyrus Harrison, Wed Nov 18 11:56:37 PST 2009
*  Declare "cmd" at top of function for strict c compile.
*
******************************************************************************/

static void
visit_add_sync(void (*cb)(void*), void *cbdata)
{
    char cmd[30];

    if(cb != NULL && callbacks->control.execute_command != NULL)
    {
        SyncCallback *sync = visit_get_sync();
        sync->id = visit_sync_id++;
        sync->cb = cb;
        sync->cbdata = cbdata;
        sprintf(cmd, "INTERNALSYNC %d", sync->id);

        (*callbacks->control.execute_command)(engine, cmd);
    }
}

/******************************************************************************
*
* Name: visit_do_sync
*
* Purpose: This function executes a callback function for the given sync id
*          and then removes the callback from the sync table.
*
* Programmer: Brad Whitlock
* Date:       Thu Mar 26 13:28:11 PDT 2009
*
* Modifications:
*
******************************************************************************/

static void
visit_do_sync(int id)
{
    SyncCallback *sync = visit_get_sync2(id);
    if(sync != NULL)
    {
        (*sync->cb)(sync->cbdata);
        sync->id = 0;
    }
}

/******************************************************************************
*
* Name: visit_handle_command_callback
*
* Purpose: This function gets installed as the engine's command callback and
*          it lets us intercept sync messages from the viewer. If we did not
*          get a sync message then we forward the command to the user-provided
*          command callback.
*
* Programmer: Brad Whitlock
* Date:       Thu Mar 26 13:28:11 PDT 2009
*
* Modifications:
*
******************************************************************************/

static void
visit_handle_command_callback(const char *cmd, const char *args, void *cbdata)
{
    if(strcmp(cmd, "INTERNALSYNC") == 0)
    {
        int id = -1;
        if(sscanf(args, "%d", &id) == 1)
            visit_do_sync(id);
    }
    else if(visit_command_callback != NULL)
    {
        (*visit_command_callback)(cmd, args, visit_command_callback_data);
    }
}

/******************************************************************************
*
* Name: visit_process_engine_command
*
* Purpose: This function processes commands from the viewer on all processors.
*          We use this function to help implement synchronization with the 
*          viewer.
*
* Programmer: Brad Whitlock
* Date:       Thu Mar 26 13:28:11 PDT 2009
*
* Modifications:
*
******************************************************************************/

#define VISIT_COMMAND_PROCESS 0
#define VISIT_COMMAND_SUCCESS 1
#define VISIT_COMMAND_FAILURE 2

static int
visit_process_engine_command(void)
{
    int command;
    LIBSIM_API_ENTER(visit_process_engine_command); 
    if(isParallel)
    {
        if (parallelRank == 0)
        {
            int success = VisItProcessEngineCommand();

            if (success)
            {
                command = VISIT_COMMAND_SUCCESS;
                BroadcastInt(&command, 0);
                return 1;
            }
            else
            {
                command = VISIT_COMMAND_FAILURE;
                BroadcastInt(&command, 0);
                return 0;
            }
        }
        else
        {
            /* Note: only through the SlaveProcessCallback callback
             * above can the rank 0 process send a VISIT_COMMAND_PROCESS
             * instruction to the non-rank 0 processes. */
            while (1)
            {
                BroadcastInt(&command, 0);
                switch (command)
                {
                case VISIT_COMMAND_PROCESS:
                    VisItProcessEngineCommand();
                    break;
                case VISIT_COMMAND_SUCCESS:
                    return 1;
                case VISIT_COMMAND_FAILURE:
                    return 0;
                }
            }
        }
    }

    command = VisItProcessEngineCommand() ? 1 : 0;
    LIBSIM_API_LEAVE(visit_process_engine_command); 
    return command;
}

/*******************************************************************************
*
* Name: CreateRandomSecurityKey
*
* Purpose: Create a random key that clients will need to connect to us.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Fri Jul 25 15:19:25 PDT 2008
*   Trace information.
*
*******************************************************************************/
static void CreateRandomSecurityKey(void)
{
    int len = 8;
    int i;
    securityKey[0] = '\0';

    LIBSIM_API_ENTER(CreateRandomSecurityKey);

#if defined(_WIN32)
    srand((unsigned)time(0));
#else
    srand48((long)(time(0)));
#endif
    for (i=0; i<len; i++)
    {
        char str[3];
#if defined(_WIN32)
        double d = (double)(rand()) / (double)(RAND_MAX);
        sprintf(str, "%02x", (int)(d * 255.));
#else
        sprintf(str, "%02x", (int)(lrand48() % 256));
#endif
        strcat(securityKey, str);
    }

    LIBSIM_API_LEAVE1(CreateRandomSecurityKey, "securityKey=%s", securityKey);
}

/*******************************************************************************
*
* Name: ReceiveSingleLineFromSocket
*
* Purpose: Receive a single line character transmission from the socket.
*          Note that this assumes it is not part of a larger transmission.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Fri Jul 25 11:54:12 PDT 2008
*   Changed some assignments to NULL to remove warnings.
*
*******************************************************************************/
static void ReceiveSingleLineFromSocket(char *buffer, size_t maxlen, VISIT_SOCKET desc)
{
    char *buf = buffer;
    char *ptr = buffer;
    char *tmp = NULL;
    int n;

    LIBSIM_API_ENTER2(ReceiveSingleLineFromSocket,
                      "maxlen=%d, desc=%d", 
                      (int)maxlen,desc);

    strcpy(buffer, "");
    tmp = strstr(buf, "\n");
    while (!tmp)
    {
        n = recv(desc, (void*)ptr, maxlen, 0);
        ptr += n;
        *ptr = '\0';
        tmp = strstr(buf, "\n");
    }
    *tmp = '\0';

    LIBSIM_API_LEAVE1(ReceiveSingleLineFromSocket,
                     "buffer=%s", buffer);
}

/*******************************************************************************
*
* Name: ReceiveContinuousLineFromSocket
*
* Purpose: Receive a single line as part of a larger transmission.  Note that
*          it assumes buffer is initialized for the first call to the empty
*          string, and that it retains its values as an intermediate buffer
*          between calls.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Fri Jul 25 11:55:05 PDT 2008
*   Changed some assignments to NULL to remove warnings.
*   
*******************************************************************************/
static char *ReceiveContinuousLineFromSocket(char *buffer, size_t maxlen, VISIT_SOCKET desc)
{
    char *buf = buffer;
    char *ptr = buffer;
    char *tmp = NULL;
    int n;

    LIBSIM_API_ENTER2(ReceiveContinuousLineFromSocket,
                      "maxlen=%d, desc=%d", 
                      (int)maxlen,desc);

    tmp = strstr(buf, "\n");
    while (!tmp)
    {
        n = recv(desc, (void*)ptr, maxlen, 0);
        ptr += n;
        *ptr = '\0';
        tmp = strstr(buf, "\n");
    }
    *tmp = '\0';

    LIBSIM_API_LEAVE1(ReceiveContinuousLineFromSocket,
                     "return %s", buffer);
    return tmp+1;
}

/*******************************************************************************
*
* Name: SendStringOverSocket
*
* Purpose: Send a single null-terminated string over a socket.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 15:14:54 PDT 2008
*  Trace information.
*
*******************************************************************************/
static int SendStringOverSocket(char *buffer, VISIT_SOCKET desc)
{
    size_t      nleft, nwritten;
    const char *sptr;
    LIBSIM_API_ENTER2(SendStringOverSocket, "buffer=%s, desc=%d", buffer,desc);
    /* Send it! */
    sptr = (const char*)buffer;
    nleft = strlen(buffer);
    while (nleft >= 1)
    {
        if ((nwritten = send(desc, (const char *)sptr, nleft, 0)) == 0)
        {
            LIBSIM_API_LEAVE1(SendStringOverSocket, 
                             "send() returned 0. return %d", TRUE);
            return FALSE;
        }
        nleft -= nwritten;
        sptr  += nwritten;
    }
    LIBSIM_API_LEAVE1(SendStringOverSocket, "return %d", TRUE);

    return TRUE;
}

/*******************************************************************************
*
* Name: BroadcastInt
*
* Purpose: Call the broadcast int callback function to broadcast an int among
*          processors.
*
* Author: Brad Whitlock, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*
*******************************************************************************/

static int
BroadcastInt(int *value, int sender)
{
    int retval = 0;

    if(sender==0)
    {
        LIBSIM_API_ENTER2(BroadcastInt, "value=%d, sender=%d", *value, sender);
    }
    else
    {
        LIBSIM_API_ENTER1(BroadcastInt, "sender=%d", sender);
    }

    if(BroadcastInt_internal != NULL)
        retval = (*BroadcastInt_internal)(value, sender);
    else
    {
        LIBSIM_MESSAGE("BroadcastInt function not set.");
    }

    LIBSIM_API_LEAVE1(BroadcastInt, "return %d", retval);
    return retval;
}

/*******************************************************************************
*
* Name: BroadcastString
*
* Purpose: Call the broadcast string callback function to broadcast a string
*          among processors.
*
* Author: Brad Whitlock, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*
*******************************************************************************/
static int
BroadcastString(char *str, int len, int sender)
{
    int retval = 0;

    if(sender==0)
    {
        LIBSIM_API_ENTER3(BroadcastString, "str=%s, len=%d, sender=%d",
                          str, len, sender);
    }
    else
    {
        LIBSIM_API_ENTER2(BroadcastString, "len=%d, sender=%d", len, sender);
    }

    if(BroadcastString_internal != NULL)
        retval = (*BroadcastString_internal)(str, len, sender);
    else
    {
        LIBSIM_MESSAGE("BroadcastString function not set.");
    }

    LIBSIM_API_LEAVE1(BroadcastString, "return %d", retval);
    return retval;
}

/*******************************************************************************
*
* Name: VerifySecurityKeys
*
* Purpose: Receive a security key over the socket and compare it to ours,
*          sending the result of the comparison back to the client.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 11:57:40 PDT 2008
*  Added casts to eliminate warnings. Added trace information.
*
*******************************************************************************/
static int VerifySecurityKeys(VISIT_SOCKET desc)
{
   int securityKeyLen;
   int offeredKeyLen;
   char offeredKey[2000] = "";

   LIBSIM_API_ENTER1(VerifySecurityKeys, "desc=%d", desc);

   if (parallelRank == 0)
   {
      /* The first thing the VCL sends is the key */
      ReceiveSingleLineFromSocket(offeredKey, 2000, desc);

      if (isParallel)
      {
         /* Broadcast the known key */
         securityKeyLen = (int)strlen(securityKey);
         BroadcastInt(&securityKeyLen, 0);
         BroadcastString(securityKey,  securityKeyLen+1, 0);

         /* Broadcast the received key */
         offeredKeyLen = (int)strlen(offeredKey);
         BroadcastInt(&offeredKeyLen, 0);
         BroadcastString(offeredKey,  offeredKeyLen+1, 0);
      }

      /* Make sure the keys match, and send a response */
      if (strcmp(securityKey, offeredKey) != 0)
      {
         SendStringOverSocket("failure\n", desc);
         LIBSIM_API_LEAVE1(VerifySecurityKeys, "return %d", FALSE);
         return FALSE;
      }
      else
      {
         SendStringOverSocket("success\n", desc);
      }
   }
   else
   {
      /* Receive the security keys and make sure they match */
      BroadcastInt(&securityKeyLen, 0);
      BroadcastString(securityKey, securityKeyLen+1, 0);
      BroadcastInt(&offeredKeyLen, 0);
      BroadcastString(offeredKey, offeredKeyLen+1, 0);

      if (strcmp(securityKey, offeredKey) != 0)
      {
         /* Error: keys didn't match */
         LIBSIM_API_LEAVE1(VerifySecurityKeys, "return %d", FALSE);
         return FALSE;
      }
   }

   LIBSIM_API_LEAVE1(VerifySecurityKeys, "return %d", TRUE);
   return TRUE;
}

/*******************************************************************************
*
* Name: GetConnectionParameters
*
* Purpose: Receive the command line to be used to initialize the VisIt engine.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Fri Jul 25 15:11:11 PDT 2008
*   Trace information.
*
*******************************************************************************/
static int GetConnectionParameters(VISIT_SOCKET desc)
{
   char buf[2000] = "";
   char *tmpbuf;
   char *nxtbuf;
   int i;

   LIBSIM_API_ENTER1(GetConnectionParameters, "desc=%d", desc);

   engine_argv = (char**)malloc(sizeof(char*) * 100);

   if (parallelRank == 0)
   {
      /* Receive the ARGV over the socket */
      engine_argc = 0;

      tmpbuf = buf;
      while (1)
      {
         nxtbuf = ReceiveContinuousLineFromSocket(tmpbuf, 2000, desc);

         if (strlen(tmpbuf) == 0)
            break;

         engine_argv[engine_argc] = strdup(tmpbuf);
         engine_argc++;
         tmpbuf = nxtbuf;
      }

      /* Broadcast them to the other processors if needed */
      if (isParallel)
      {
         BroadcastInt(&engine_argc, 0);
         for (i = 0 ; i < engine_argc; i++)
         {
            int len = strlen(engine_argv[i]);
            BroadcastInt(&len, 0);
            BroadcastString(engine_argv[i], len+1, 0);
         }
      }
   }
   else
   {
      /* Receive the ARGV */
      BroadcastInt(&engine_argc, 0);
      for (i = 0 ; i < engine_argc; i++)
      {
         int len;
         BroadcastInt(&len, 0);
         BroadcastString(buf, len+1, 0);
         engine_argv[i] = strdup(buf);
      }
   }

   LIBSIM_API_LEAVE1(GetConnectionParameters, "return %d", TRUE);

   return TRUE;
}

/*******************************************************************************
*
* Name: CreateEngineAndConnectToViewer
*
* Purpose: Intialize the engine and connect to the viewer.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Fri Jul 25 14:33:37 PDT 2008
*   Trace information.
*
*******************************************************************************/

static int CreateEngineAndConnectToViewer(void)
{
    LIBSIM_API_ENTER(CreateEngineAndConnectToViewer);

    /* get the engine */
    LIBSIM_MESSAGE("Calling visit_engine");
    engine = (*callbacks->control.get_engine)();
    if (!engine)
    {
        LIBSIM_API_LEAVE1(CreateEngineAndConnectToViewer,
                         "engine could not be allocated. return %d",
                         FALSE);
        return FALSE;
    }

    LIBSIM_MESSAGE_STRINGLIST("Calling visit_initialize: argv",
                              engine_argc, engine_argv);
    if (!(*callbacks->control.initialize)(engine, engine_argc, engine_argv))
    {
        VisItDisconnect();
        LIBSIM_API_LEAVE1(CreateEngineAndConnectToViewer,
                         "visit_initialize failed. return %d",
                         FALSE);
        return FALSE;
    }

    LIBSIM_MESSAGE_STRINGLIST("Calling visit_connectviewer: argv",
                              engine_argc, engine_argv);
    if (!(*callbacks->control.connect_viewer)(engine, engine_argc, engine_argv))
    {
        VisItDisconnect();
        LIBSIM_API_LEAVE1(CreateEngineAndConnectToViewer,
                         "visit_connectviewer failed. return %d", 
                         FALSE);
        return FALSE;
    }

    LIBSIM_API_LEAVE1(CreateEngineAndConnectToViewer,"return %d", TRUE);
    return TRUE;
}

/*******************************************************************************
*
* Name: GetLocalhostName
*
* Purpose: Determine the true local host name.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 14:41:30 PDT 2008
*  Trace information.
*
*  Brad Whitlock, Fri Nov 12 23:40:23 PST 2010
*  I added an extra layer of gethostbyname to match RemoteProcess.
*
*******************************************************************************/
static int GetLocalhostName(void)
{
    char localhostStr[256];
    struct hostent *localhostEnt = NULL;

    LIBSIM_API_ENTER(GetLocalhostName);

    LIBSIM_MESSAGE("Calling gethostname");
    if (gethostname(localhostStr, 256) == -1)
    {
        /* Couldn't get the hostname, it's probably invalid */
        LIBSIM_API_LEAVE1(GetLocalhostName, 
                         "gethostname failed. return=%d", FALSE);
        return FALSE;
    }
    LIBSIM_MESSAGE1("gethostname returned %s\n", localhostStr);

    LIBSIM_MESSAGE("Calling gethostbyname");
    localhostEnt = gethostbyname(localhostStr);
    if (localhostEnt == NULL)
    {
        /* Couldn't get the full host entry; it's probably invalid */
        LIBSIM_MESSAGE("gethostbyname failed. call gethostbyname(localhost)");
        
        strcpy(localhostStr, "localhost");
        localhostEnt = gethostbyname(localhostStr);
        if(localhostEnt != NULL)
        {
            sprintf(localhost, "%s", localhostEnt->h_name);
            LIBSIM_MESSAGE1("Using return value of gethostbyname: %s", localhost);
        }
        else
        {
            LIBSIM_MESSAGE("Use \"localhost\" as the host name.");
            strcpy(localhost, "localhost");
        }
    }
    else
    {
        sprintf(localhost, "%s", localhostEnt->h_name);
        LIBSIM_MESSAGE1("Using return value of gethostbyname: %s", localhost);
    }
    LIBSIM_API_LEAVE1(GetLocalhostName, "return %s", localhost);

    return TRUE;
}

/*******************************************************************************
*
* Name: StartListening
*
* Purpose: Find a port and start listening.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*
*******************************************************************************/
static int StartListening(void)
{
    int portFound = FALSE;
    int on = 1;
    int err;

    LIBSIM_API_ENTER(StartListening);
    listenSocket = socket(AF_INET, SOCK_STREAM, 0);
    if (listenSocket < 0)
    {
        /* Cannot open a socket. */
        LIBSIM_API_LEAVE1(StartListening, "socket() failed. return=%d", FALSE);
        return FALSE;
    }

    /*
     * Look for a port that can be used.
     */
    LIBSIM_MESSAGE("Looking for a listening port");
    listenSockAddr.sin_family = AF_INET;
    listenSockAddr.sin_addr.s_addr = htonl(INADDR_ANY);
    listenPort = INITIAL_PORT_NUMBER;
    while (!portFound && listenPort < 32767)
    {
        listenSockAddr.sin_port = htons(listenPort);
#if !defined(_WIN32)
        setsockopt(listenSocket, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on));
#endif

        err = bind(listenSocket, (struct sockaddr *)&listenSockAddr,
                   sizeof(listenSockAddr));
        if (err)
        {
            listenPort++;
        }
        else
        {
            portFound = TRUE;
        }
    }
    if (!portFound)
    {
        /* Cannot find unused port. */
       listenSocket = VISIT_INVALID_SOCKET;
       LIBSIM_API_LEAVE1(StartListening, "port not found. return %d", FALSE);
       return FALSE;
    }

    LIBSIM_MESSAGE("Calling listen() on socket");
    err = listen(listenSocket, 5);
    if (err)
    {
       listenSocket = VISIT_INVALID_SOCKET;
       LIBSIM_API_LEAVE1(StartListening, "listen() failed. return %d", FALSE);
       return FALSE;
    }

    LIBSIM_API_LEAVE1(StartListening, "return %d", TRUE);
    return TRUE;
}

/*******************************************************************************
*
* Name: AcceptConnection
*
* Purpose: Perform the accept on the listen socket.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 14:47:55 PDT 2008
*  Trace information.
*
*******************************************************************************/
static VISIT_SOCKET AcceptConnection(void)
{
    VISIT_SOCKET desc = VISIT_INVALID_SOCKET;
    int opt = 1;

    LIBSIM_API_ENTER(AcceptConnection);

    /* Wait for the socket to become available on the other side. */
    do
    {
#ifdef HAVE_SOCKLEN_T
        socklen_t len;
#else
#ifdef __APPLE__
        unsigned int len;
#else
        int len;
#endif
#endif
        len = sizeof(struct sockaddr);
        LIBSIM_MESSAGE("Calling accept()");
        desc = accept(listenSocket, (struct sockaddr *)&listenSockAddr, &len);
    }
    while (desc == VISIT_INVALID_SOCKET
#ifndef _WIN32
           && errno != EWOULDBLOCK
#endif
           );

    /* Disable Nagle algorithm. */
#if defined(_WIN32)
    setsockopt(desc, IPPROTO_TCP, TCP_NODELAY,
               (const char FAR*)&opt, sizeof(int));
#else
    setsockopt(desc, IPPROTO_TCP, TCP_NODELAY, &opt, sizeof(int));
#endif
    
    LIBSIM_API_LEAVE1(AcceptConnection, "desc=%d", desc);
    return desc;
}

/*******************************************************************************
*
* Name: GetHomeDirectory
*
* Purpose: Return the true home directory path.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*
*******************************************************************************/
static const char *GetHomeDirectory(void)
{
#ifdef _WIN32
    return "C:/Users/Brad";
#else
    struct passwd *users_passwd_entry = NULL;

    LIBSIM_API_ENTER(GetHomeDirectory);
    users_passwd_entry = getpwuid(getuid());
    LIBSIM_API_LEAVE1(GetHomeDirectory, "homedir=%s", 
                     users_passwd_entry->pw_dir);

    return users_passwd_entry->pw_dir;
#endif
}

/*******************************************************************************
*
* Name: EnsureSimulationDirectoryExists
*
* Purpose: Make the simulations directory.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Fri Jul 25 14:56:47 PDT 2008
*   Trace information.
*
*******************************************************************************/
static void EnsureSimulationDirectoryExists(void)
{
    char str[1024];
    LIBSIM_API_ENTER(EnsureSimulationDirectoryExists);

    SNPRINTF(str, 1024, "%s/.visit", GetHomeDirectory());
    VisItMkdir(str, 7*64 + 7*8 + 7);
    LIBSIM_MESSAGE1("mkdir %s", str);

    SNPRINTF(str, 1024, "%s/.visit/simulations", GetHomeDirectory());
    VisItMkdir(str, 7*64 + 7*8 + 7);
    LIBSIM_MESSAGE1("mkdir %s", str);

    LIBSIM_API_LEAVE(EnsureSimulationDirectoryExists);
}

/*******************************************************************************
*
* Name: RemoveSimFile
*
* Purpose: This will delete the sim file from the file system.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 14:57:51 PDT 2008
*  Trace information.
*
*******************************************************************************/
static void RemoveSimFile(void)
{
    LIBSIM_API_ENTER(RemoveSimFile);
    LIBSIM_MESSAGE1("unlink(%s)", simulationFileName);
    unlink(simulationFileName);
    LIBSIM_API_LEAVE(RemoveSimFile);
}

/*******************************************************************************
*
* Name: visit_get_runtime_function
*
* Purpose: Get the function pointer from the runtime library
*
* Author: Brad Whitlock
*
* Modifications:
*   Brad Whitlock, Thu Nov 25 23:26:23 PST 2010
*   Ported to Windows.
*
*******************************************************************************/

void *
visit_get_runtime_function(const char *name)
{
    void *f = NULL;

    LIBSIM_API_ENTER2(visit_get_runtime_function, "handle=%p, name=%s", dl_handle, name);
#ifdef _WIN32
    f = GetProcAddress(dl_handle, name);
#else
    f = dlsym(dl_handle, name);
#endif
    if(f == NULL)
    {
        sprintf(lastError, "Function %s could not be located in the runtime.\n", name);
    }
    LIBSIM_API_LEAVE1(visit_get_runtime_function, "func=%p", (void*)f);
    return f;
}

/*******************************************************************************
*
* Name: LoadVisItLibrary_XXX
*
* Purpose: Loads the Visit runtime library, setting dl_handle and lastError.
*
* Author: Brad Whitlock
*
* Modifications:
*
*******************************************************************************/

#ifdef _WIN32
static int LoadVisItLibrary_Windows(void)
{
    char visitpath[1024], lib[1024];

    LIBSIM_API_ENTER(LoadVisItLibrary_Windows);
    GetVisItDirectory(visitpath, 1024);

    /* load library */
    if (isParallel)
    {
        SNPRINTF(lib, 256, "%s\\simV2runtime_par.dll", visitpath);
    }
    else
    {
        SNPRINTF(lib, 256, "%s\\simV2runtime_ser.dll", visitpath);
    }

    LIBSIM_MESSAGE1("Loading runtime: %s", lib);
    SetErrorMode(0);
    dl_handle = LoadLibrary(lib);

    if (dl_handle == NULL)
    {
        WCHAR msg[1024];
        va_list *va = 0;
LIBSIM_MESSAGE1("Error: %d", GetLastError());
LIBSIM_MESSAGE1("Error: %x", GetLastError());
        FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM, 0, GetLastError(), 0,
                      msg, 1024, va);
        wprintf(L"Formatted message: %s\n", msg);

        /*SNPRINTF(lastError, 1024, "Failed to open the VisIt library: %s\n", msg);*/
    }

    LIBSIM_API_LEAVE(LoadVisItLibrary_Windows);

    return (dl_handle != NULL) ? VISIT_OKAY : VISIT_ERROR;
}
#else
static int LoadVisItLibrary_UNIX(void)
{
    char lib[256];

    /* load library */
#ifdef __APPLE__
    const char *extension = "dylib";
    const char *LD_LIBRARY_PATH = "DYLD_LIBRARY_PATH";
#else
    const char *extension = "so";
    const char *LD_LIBRARY_PATH = "LD_LIBRARY_PATH";
#endif

    if (isParallel)
    {
        sprintf(lib, "libsimV2runtime_par.%s", extension);
    }
    else
    {
        sprintf(lib, "libsimV2runtime_ser.%s", extension);
    }

    LIBSIM_MESSAGE1("Calling dlopen(%s)", lib);
    dl_handle = dlopen(lib, RTLD_NOW | RTLD_GLOBAL);

    if (dl_handle == NULL)
    {
        LIBSIM_MESSAGE1("Failed to open the VisIt library: %s\n", dlerror())
        LIBSIM_MESSAGE1("Calling getenv(%s)", LD_LIBRARY_PATH);

        if(getenv(LD_LIBRARY_PATH) != NULL)
        {
            char *libpath = NULL, *start = NULL, *ptr = NULL;
            libpath = start = strdup(getenv(LD_LIBRARY_PATH));
            LIBSIM_MESSAGE1("getenv returned: %s", getenv(LD_LIBRARY_PATH));

            /* Iterate over all paths in the LD_LIBRARY_PATH and try each one
             * until we get a library to open.
             */
            do
            {
                /* If there is a separator then null it out*/
                ptr = strstr(libpath, ":");
                if (ptr)
                    *ptr = '\0';

                /* Try to use the current libpath to open the library */
                if (isParallel) 
                    sprintf(lib, "%s/libsimV2runtime_par.%s", libpath, extension);
                else
                    sprintf(lib, "%s/libsimV2runtime_ser.%s", libpath, extension);
                LIBSIM_MESSAGE1("Calling dlopen(%s)", lib);
                dl_handle = dlopen(lib, RTLD_NOW | RTLD_GLOBAL);
                if(dl_handle == NULL)
                {
                     LIBSIM_MESSAGE1("dlopen error: %s", dlerror());
                }
                if(ptr != NULL)
                    libpath = ptr + 1;
                else
                    libpath += strlen(libpath);
            } while(dl_handle == NULL && *libpath != '\0');
            free(start);
        }
    }

    if (dl_handle == NULL)
        sprintf(lastError, "Failed to open the VisIt library: %s\n", dlerror());
    
    return (dl_handle != NULL) ? VISIT_OKAY : VISIT_ERROR;
}
#endif

/*******************************************************************************
*
* Name: CloseVisItLibrary
*
* Purpose: Closes the Visit runtime library, setting dl_handle.
*
* Author: Brad Whitlock
*
* Modifications:
*
*******************************************************************************/

static void CloseVisItLibrary(void)
{
#ifdef _WIN32
    if(dl_handle != NULL)
    {
        FreeLibrary(dl_handle);
        dl_handle = NULL;
    }
#else
    /* Call dlclose(dl_handle) ???*/
    dl_handle = NULL;
#endif
}

/*******************************************************************************
*
* Name: LoadVisItLibrary
*
* Purpose: Load the DLL/SO of the VisIt Engine and get the needed function
*          pointers from it.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Thu Jan 25 14:58:26 PST 2007
*   Added update plots and execute_command.
*
*   Jeremy Meredith, Fri Nov  2 18:06:42 EDT 2007
*   Use dylib as the extension for OSX.
*
*   Brad Whitlock, Fri Feb 13 09:54:47 PST 2009
*   I added code to get data callback setting functions from the library.
*
*   Brad Whitlock, Tue Feb 16 15:09:52 PST 2010
*   Try a list of libraries in the library path.
*
*   Brad Whitlock, Thu Nov 25 23:03:23 PST 2010
*   I moved the code to load libraries into helper functions to separate out
*   some Windows-only logic.
*
*******************************************************************************/

static int LoadVisItLibrary(void)
{
    int status = VISIT_ERROR;
    LIBSIM_API_ENTER(LoadVisItLibrary);

    /* Load the library */
#ifdef _WIN32
    status = LoadVisItLibrary_Windows();
#else
    status = LoadVisItLibrary_UNIX();
#endif
    if (status == VISIT_ERROR)
    {
        LIBSIM_API_LEAVE2(LoadVisItLibrary, "%s: return %d", lastError, FALSE);
        return FALSE;
    }

    callbacks = (visit_callback_t *)malloc(sizeof(visit_callback_t));
    memset(callbacks, 0, sizeof(visit_callback_t));

#define SAFE_DLSYM(D,F,T,N,DECORATE) \
    callbacks->D.F = (DECORATE(T))visit_get_runtime_function(N); \
    if (!callbacks->D.F) \
    { \
        CloseVisItLibrary();\
        LIBSIM_API_LEAVE2(LoadVisItLibrary, "%s: return %d", lastError, FALSE); \
        return FALSE; \
    }
#define NO_DECORATE(A) A
#define SET_DECORATE(A) void (*)(A,void*)
#define QUOTED(A) #A
#define CONTROL_DLSYM(F,T) SAFE_DLSYM(control,F,T,QUOTED(simv2_##F),NO_DECORATE)
#define DATA_DLSYM(F,T)    SAFE_DLSYM(data,set_##F,T,QUOTED(simv2_set_##F),SET_DECORATE)

   /* Get the control functions fom the library. */
   CONTROL_DLSYM(get_engine,                 void *(*)(void));
   CONTROL_DLSYM(get_descriptor,             int (*)(void *));
   CONTROL_DLSYM(process_input,              int (*)(void *));
   CONTROL_DLSYM(initialize,                 int (*)(void *, int, char **));
   CONTROL_DLSYM(connect_viewer,             int (*)(void *, int, char **));
   CONTROL_DLSYM(time_step_changed,          void (*)(void *));
   CONTROL_DLSYM(execute_command,            void (*)(void *,const char*));
   CONTROL_DLSYM(disconnect,                 void (*)());
   CONTROL_DLSYM(set_slave_process_callback, void (*)(void (*)()));
   CONTROL_DLSYM(set_command_callback,       void (*)(void*,void (*)(const char*,const char*,void*),void*));
   CONTROL_DLSYM(save_window,                int (*)(void*,const char *,int,int,int));
   CONTROL_DLSYM(debug_logs,                 void (*)(int,const char *));

   /* Get the data functions from the library. */
   DECLARE_DATA_CALLBACKS(DATA_DLSYM, NO_ARGS)

#ifdef VISIT_QUERIES_VIA_LIBSIM
   /* Load the VisIt query functions.*/
   LoadVisItQueries();
#endif

   LIBSIM_API_LEAVE1(LoadVisItLibrary, "return %d", TRUE);
   return TRUE;
}

/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
                             Public API Functions
 *******************************************************************************
 *******************************************************************************
 ******************************************************************************/


/*******************************************************************************
*
* Name: VisItSetBroadcastIntFunction
*
* Purpose: Set the callback to broadcast an integer.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Fri Jul 25 15:33:27 PDT 2008
*   Trace information.
*
*******************************************************************************/
void  VisItSetBroadcastIntFunction(int (*bi)(int *, int))
{
   LIBSIM_API_ENTER1(VisItSetBroadcastIntFunction, "func=%p", (void*)bi);
   BroadcastInt_internal = bi;
   LIBSIM_API_LEAVE(VisItSetBroadcastIntFunction);
}

/*******************************************************************************
*
* Name: VisItSetBroadcastStringFunction
*
* Purpose: Set the callback to broadcast a string.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 15:33:37 PDT 2008
*  Trace information.
*
*******************************************************************************/
void  VisItSetBroadcastStringFunction(int (*bs)(char *, int, int))
{
   LIBSIM_API_ENTER1(VisItSetBroadcastStringFunction, "func=%p", (void*)bs);
   BroadcastString_internal = bs;
   LIBSIM_API_LEAVE(VisItSetBroadcastStringFunction);
}

/*******************************************************************************
*
* Name: VisItSetParallel
*
* Purpose: Set whether or not we have to work in parallel.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 15:33:37 PDT 2008
*  Trace information.
*
*******************************************************************************/
void  VisItSetParallel(int ispar)
{
   LIBSIM_API_ENTER1(VisItSetParallel, "isparallel=%d", ispar);
   isParallel = ispar;
   LIBSIM_API_LEAVE(VisItSetParallel);
}

/*******************************************************************************
*
* Name: VisItSetParallelRank
*
* Purpose: Set the rank of this processor.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 15:33:37 PDT 2008
*  Trace information.
*
*******************************************************************************/
void  VisItSetParallelRank(int pr)
{
   LIBSIM_API_ENTER1(VisItSetParallelRank, "rank=%d", pr);
   parallelRank = pr;
   LIBSIM_API_LEAVE(VisItSetParallelRank);
}

/*******************************************************************************
*
* Name: VisItSetDirectory
*
* Purpose: Set the top-level directory for VisIt.  This can either be a
*          installed or development version.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 15:33:37 PDT 2008
*  Trace information.
*
*******************************************************************************/
void VisItSetDirectory(char *d)
{
   LIBSIM_API_ENTER1(VisItSetDirectory, "dir=%s", d);
   if (!visit_directory)
      visit_directory = (char*)(malloc(1000));
   strcpy(visit_directory, d);
   LIBSIM_API_LEAVE(VisItSetDirectory);
}

/*******************************************************************************
*
* Name: VisItSetOptions
*
* Purpose: Set command-line options needed to pick up the right VisIt.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*
*******************************************************************************/
void VisItSetOptions(char *o)
{
   LIBSIM_API_ENTER1(VisItSetOptions, "options=%s", o);
   if (!visit_options)
      visit_options = (char*)(malloc(1000));
   strcpy(visit_options, o);
   LIBSIM_API_LEAVE(VisItSetOptions);
}

/*******************************************************************************
*
* Name: VisItSetupEnvironment
*
* Purpose: Try to determine the environment variables that the VisIt Engine
*          needs to run.  The VisIt script can tell us this.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Eric Brugger, Thu Sep 14 12:53:46 PDT 2006
*   Changed the routine to allocate the buffer new_env to be ENV_BUF_SIZE
*   bytes in size instead of 10000.
*
*   Brad Whitlock, Fri Jul 25 15:37:05 PDT 2008
*   Trace information.
*
*   Brad Whitlock, Thu Nov 25 23:40:34 PST 2010
*   Windows port.
*
*******************************************************************************/
int VisItSetupEnvironment(void)
{
#ifdef _WIN32
    WORD wVersionRequested;
    WSADATA wsaData;
    char visitpath[1024], tmp[1024];

    LIBSIM_API_ENTER(VisItSetupEnvironment);
    GetVisItDirectory(visitpath, 1024);

    /* Tell Windows that we want to get DLLs from this path */
    SetDllDirectory(visitpath);

    /* Set the VisIt home dir. */
    SNPRINTF(tmp, 1024, "VISITHOME=%s", visitpath);
    LIBSIM_MESSAGE(tmp);
    putenv(tmp);

    /* Set the plugin dir. */
    SNPRINTF(tmp, 1024, "VISITPLUGINDIR=%s", visitpath);
    LIBSIM_MESSAGE(tmp);
    putenv(tmp);

    // Initiate the use of a Winsock DLL (WS2_32.DLL), necessary for sockets.
    wVersionRequested = MAKEWORD(2,2);
    WSAStartup(wVersionRequested, &wsaData);
#else
   char *new_env = (char*)(malloc(ENV_BUF_SIZE));
   int done = 0;
   char *ptr;

   LIBSIM_API_ENTER(VisItSetupEnvironment);

   /* Try the one specified in by the visit_dir command first */
   if (visit_directory)
   {
      char path[200];
      sprintf(path, "%s/bin/visit", visit_directory);
      done = ReadEnvironmentFromCommand(path, new_env);
   }

   /* Try the one in their path next */
   if (!done)
   {
      done = ReadEnvironmentFromCommand("visit", new_env);
   }

   /* If we still can't find it, try the one in /usr/gapps/visit */
   if (!done)
   {
      done = ReadEnvironmentFromCommand("/usr/gapps/visit/bin/visit", new_env);
   }

   if (!done)
   {
      LIBSIM_API_LEAVE1(VisItSetupEnvironment, "return %d", FALSE);
      return FALSE;
   }

   /* Do a bunch of putenv calls; it should already be formatted correctly */
   ptr = new_env;
   while (ptr[0]!='\0')
   {
      int i = 0;
      while (ptr[i]!='\n')
         i++;
      ptr[i] = '\0';

      LIBSIM_MESSAGE1("putenv(%s)", ptr);
      putenv(ptr);

      ptr += i+1;
   }
   /* free(new_env); <--- NO!  You are not supposed to free this memory! */
#endif
   LIBSIM_API_LEAVE1(VisItSetupEnvironment, "return %d", TRUE);
   return TRUE;
}

/*******************************************************************************
*
* Name: VisItInitializeSocketAndDumpSimFile
*
* Purpose: Start listening, and write the file telling clients how to connect.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Shelly Prevost, Wed Jan 25 08:50:44 PST 2006
*   Added the guifile argument.
*
*   Brad Whitlock, Fri Jul 25 15:40:03 PDT 2008
*   Trace information.
*
*******************************************************************************/
int VisItInitializeSocketAndDumpSimFile(const char *name,
                                        const char *comment,
                                        const char *path,
                                        const char *inputfile,
                                        const char *guifile,
                                        const char *absoluteFilename)
{
    FILE *file = NULL;
    LIBSIM_API_ENTER(VisItInitializeSocketAndDumpSimFile)
    LIBSIM_MESSAGE1("name=%s", name);
    LIBSIM_MESSAGE1("comment=%s", comment);
    LIBSIM_MESSAGE1("path=%s", path);
    LIBSIM_MESSAGE1("inputfile=%s", inputfile);
    LIBSIM_MESSAGE1("guifile=%s", guifile);
    LIBSIM_MESSAGE1("absoluteFilename=%s", absoluteFilename);

    EnsureSimulationDirectoryExists();
    CreateRandomSecurityKey();
    
    if ( !absoluteFilename )
    {
        SNPRINTF(simulationFileName, 255, "%s/.visit/simulations/%012d.%s.sim2",
                 GetHomeDirectory(), (int)time(NULL), name);
    }
    else
    {
        SNPRINTF(simulationFileName, 255, "%s", absoluteFilename);
    }

    if (!GetLocalhostName())
    {
        LIBSIM_API_LEAVE1(VisItInitializeSocketAndDumpSimFile, "return %d", FALSE);
        return FALSE;
    }

    if (!StartListening())
    {
        LIBSIM_API_LEAVE1(VisItInitializeSocketAndDumpSimFile, "return %d", FALSE);
        return FALSE;
    }

    LIBSIM_MESSAGE1("Opening sim file %s", simulationFileName);  
    file = fopen(simulationFileName, "wt");
    if (!file)
    {
        LIBSIM_API_LEAVE1(VisItInitializeSocketAndDumpSimFile, "return %d", FALSE);
        return FALSE;
    }

    atexit(RemoveSimFile);

    fprintf(file, "host %s\n", localhost);
    fprintf(file, "port %d\n", listenPort);
    fprintf(file, "key %s\n", securityKey);
    if (path)
        fprintf(file, "path %s\n", path);
    if (inputfile)
        fprintf(file, "inputfile %s\n", inputfile);
    if (comment)
        fprintf(file, "comment %s\n", comment);
     if ( guifile )
       fprintf(file, "uiFile %s\n", guifile);

    fclose(file);

    LIBSIM_API_LEAVE1(VisItInitializeSocketAndDumpSimFile, "return %d", TRUE);
    return TRUE;
}

/*******************************************************************************
*
* Name: VisItDetectInput
*
* Purpose: Determine who needs attention: the listen socket, the client
*          socket, or the console.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Wed Jun 25 16:28:02 PDT 2008
*   Change 500000us (1/2 second) timeout to zero. The 500000 must have been
*   left in during an old debugging exercise.
*
*******************************************************************************/

int
VisItDetectInput(int blocking, int consoleFileDescriptor)
{
    return VisItDetectInputWithTimeout(blocking, 0, consoleFileDescriptor);
}

#ifdef _WIN32
static int inputThreadsStarted = 0;
static HANDLE stdinevent = 0;
static WSAEVENT listenevent = 0;
static HANDLE engineevent = 0;
static HANDLE listeneventCB = 0;
static HANDLE engineeventCB = 0;

DWORD WINAPI
stdin_thread(LPVOID param)
{
    while(1)
    {
        WaitForSingleObject(GetStdHandle(STD_INPUT_HANDLE), INFINITE);
        fprintf(stderr, "stdin thread signaling input\n");
        SetEvent(stdinevent);
    }

    return 0;
}

#define SELECT_ENGINE_SOCKET

static void
VLogWindowsSocketError()
{
    switch(WSAGetLastError())
    {
    case WSANOTINITIALISED:
        fprintf(stderr,"WSAENOTINITIALISED: WSAStartup() must be called before using this API.");
        break;
    case WSAENETDOWN:
        fprintf(stderr,"WSAENETDOWN: The network subsystem or the associated service provider has failed.");
        break;
    case WSAEAFNOSUPPORT:
        fprintf(stderr,"WSAEAFNOSUPPORT: The specified address family is not supported.");
        break;
    case WSAEINPROGRESS:
        fprintf(stderr,"WSAEINPROGRESS: A blocking Winsock 1.1 call is in progress, or the service provider is still processing a callback function.");
        break;
    case WSAEFAULT:
        fprintf(stderr,"WSAEFAULT:");
        break;
    case WSAEINTR:
        fprintf(stderr,"WSAEINTR: A blocking WinSock 1.1 call was canceled via WSACancelBlockingCall.");
        break;
    case WSAEMFILE:
        fprintf(stderr,"WSAEMFILE: No more socket descriptors are available.");
        break;
    case WSAENOBUFS:
        fprintf(stderr,"WSAENOBUFS: No buffer space is available. The socket cannot be created.");
        break;
    case WSAEPROTONOSUPPORT:
        fprintf(stderr,"WSAEPROTONOSUPPORT: The specified protocol is not supported.");
        break;
    case WSAEPROTOTYPE:
        fprintf(stderr,"WSAEPROTOTYPE: The specified protocol is the wrong type for this socket.");
        break;
    case WSAESOCKTNOSUPPORT:
        fprintf(stderr,"WSAESOCKTNOSUPPORT: The specified socket type is not supported in this address family.");
        break;
    case WSAHOST_NOT_FOUND:
        fprintf(stderr,"WSAHOST_NOT_FOUND: Authoratiative Answer Host not found.");
        break;
    case WSATRY_AGAIN:
        fprintf(stderr,"WSATRY_AGAIN: Non-Authoratative Host not found, or server failure.");
        break;
    case WSANO_RECOVERY:
        fprintf(stderr,"WSANO_RECOVERY: Non-recoverable error occurred.");
        break;
    case WSANO_DATA:
        fprintf(stderr,"WSANO_DATA: Valid name, no data record of requested type.");
        break;
    case WSAEINVAL:
        fprintf(stderr,"WSAEINVAL: See documentation for: ");
        break;
    case WSAENETRESET:
        fprintf(stderr,"WSAENETRESET: The connection has been broken due to keep-alive activity detecting a failure while the operation was in progress.");
        break;
    case WSAENOPROTOOPT:
        fprintf(stderr,"WSAENOPROTOOPT: The option is unknown or unsupported for the specified provider or socket.");
        break;
    case WSAENOTCONN:
        fprintf(stderr,"WSAENOTCONN: Connection has been reset when SO_KEEPALIVE is set.");
        break;
    case WSAENOTSOCK:
        fprintf(stderr,"WSAENOTSOCK: The descriptor is not a socket.");
        break;
    case WSAEADDRINUSE:
        fprintf(stderr,"WSAEADDRINUSE: The socket's local address space is already in use and the socket was not marked to allow address reuse.");
        break;
    case WSAEALREADY:
        fprintf(stderr,"WSAEALREADY: A non-blocking connect() call is in progress or the service provider is still processing a callback function.");
        break;
    case WSAEADDRNOTAVAIL:
        fprintf(stderr,"WSAADDRNOTAVAIL: The remote address is not valid (e.g. ADDR_ANY).");
        break;
    case WSAECONNREFUSED:
        fprintf(stderr,"WSAECONNREFUSED: The attempt to connect was forcefully rejected.");
        break;
    case WSAEISCONN:
        fprintf(stderr,"WSAEISCONN: The socket is already connected.");
        break;
    case WSAENETUNREACH:
        fprintf(stderr,"WSAENETUNREACH: The network can't be reached from this host at this time.");
        break;
    case WSAETIMEDOUT:
        fprintf(stderr,"WSAETIMEDOUT: Attempt to connect timed out without establishing a connection.");
        break;
    case WSAEWOULDBLOCK:
        fprintf(stderr,"WSAEWOULDBLOCK: The socket it marked as non-blocking and the connection cannot be completed immediately.");
        break;
    case WSAEACCES:
        fprintf(stderr,"WSAEACCES: Attempt to connect datagram socket to broadcast address failed.");
        break;
    case WSAEOPNOTSUPP:
        fprintf(stderr,"WSAENOTSUPP: The referenced socket is not of a type that supports listen().");
        break;
    default:
        fprintf(stderr, "WSA error: %d", WSAGetLastError());
    }
    fprintf(stderr, "\n");
}

DWORD WINAPI
select_thread(LPVOID param)
{
    while(engineSocket != VISIT_INVALID_SOCKET ||
          listenSocket != VISIT_INVALID_SOCKET)
    {
        fd_set readSet;
        int    ignored = 0, status = 0;
        struct timeval SmallTimeout = {0, 50000};

        FD_ZERO(&readSet);

        /* If we're connected, select on the control socket */
#ifdef SELECT_ENGINE_SOCKET
        if(engineSocket != VISIT_INVALID_SOCKET)
        {
            FD_SET(engineSocket, &readSet);
            fprintf(stderr, "select_thread: selecting engine control socket\n");
        }
#endif

        /* If we're connected, do *not* select on the listen socket */
        /* This forces us to have only one client at a time. */
        if (engineSocket == VISIT_INVALID_SOCKET &&
            listenSocket != VISIT_INVALID_SOCKET)
        {
            FD_SET(listenSocket, &readSet);
            fprintf(stderr, "select_thread: selecting listen socket for inbound connections\n");
        }

        fprintf(stderr, "select_thread: calling select\n");
        status = select(ignored, &readSet, (fd_set*)NULL, 
                        (fd_set*)NULL, NULL);

//haveEngineSocket ? &SmallTimeout : NULL);

        if(status == SOCKET_ERROR)
        {
            fprintf(stderr, "SOCKET_ERROR\n");
            VLogWindowsSocketError();
        }
        else if(status == 0)
        {
            // timed out
//fprintf(stderr, "timed out!\n");
        }
        else if(status > 0)
        {
             if (listenSocket != VISIT_INVALID_SOCKET &&
                 FD_ISSET(listenSocket, &readSet))
             {
//               fprintf(stderr, "Send a listen event!\n");
                 SetEvent(listenevent);
                 // wait for it to be done.
                 WaitForSingleObject(listeneventCB, INFINITE);
             }
#ifdef SELECT_ENGINE_SOCKET
             else if (engineSocket != VISIT_INVALID_SOCKET &&
                      FD_ISSET(engineSocket, &readSet))
             {
//               fprintf(stderr, "Send an engine event!\n");
                 SetEvent(engineevent);

                 // wait for it to be done.
                 WaitForSingleObject(engineeventCB, INFINITE);
             }
#endif
        }
        else
        {
            // error
        }
    }

    return 0;
}

int
VisItDetectInputWithTimeout(int blocking, int timeoutVal,
    int consoleFileDescriptor)
{
    HANDLE handles[3];
    int n, msec;

    LIBSIM_API_ENTER2(VisItDetectInput, "blocking=%d, consoleFile=%d",
                      blocking, consoleFileDescriptor);

    msec = timeoutVal / 1000;

    if(!inputThreadsStarted)
    {
        DWORD stdin_threadid, select_threadid;

        stdinevent  = CreateEvent(NULL, FALSE, FALSE, NULL);
        listenevent = CreateEvent(NULL, FALSE, FALSE, NULL);
        engineevent = CreateEvent(NULL, FALSE, FALSE, NULL);
        listeneventCB = CreateEvent(NULL, FALSE, FALSE, NULL);
        engineeventCB = CreateEvent(NULL, FALSE, FALSE, NULL);

#if 0
        fprintf(stderr, "Creating stdin thread\n");
        if(!CreateThread(NULL, 0, stdin_thread,
                         NULL, 0, &stdin_threadid))
        {
            LIBSIM_API_LEAVE1(VisItDetectInput,
                              "Unable to create console input thread. return %d",
                              -4);
            return -4;
        }
#endif

        fprintf(stderr, "Creating socket input thread\n");
        if(!CreateThread(NULL, 0, select_thread,
                         NULL, 0, &select_threadid))
        {
            LIBSIM_API_LEAVE1(VisItDetectInput,
                              "Unable to create socket input thread. return %d",
                              -4);
            return -4;
        }

        inputThreadsStarted = 1;
    }

    // Wait for any of these events to occur.
    handles[0] = listenevent;
    handles[1] = engineevent;
    //handles[2] = GetStdHandle(STD_INPUT_HANDLE);
    n = WaitForMultipleObjects(2, handles, FALSE, blocking ? INFINITE : msec);

    //fprintf(stderr, "MsgWaitForMultipleObjects: returned %d\n", n);

    if(n == WAIT_TIMEOUT)
    {
        LIBSIM_API_LEAVE1(VisItDetectInput,
                          "Okay - Timed out. return %d",
                          0);
        return 0;
    }
    else if(n == WAIT_OBJECT_0)
    {
        // we received an event on the listen socket. If it's an event
        // indicating that there's input then tell the user it's okay
        // to read from the listen socket.
        LIBSIM_API_LEAVE1(VisItDetectInput,
                          "WAIT_OBJECT_0: Listen  socket input. return %d",
                          1);
//        SetEvent(listeneventCB);
        return 1;
    }
    else if(n == WAIT_OBJECT_0+1)
    {
        // we received an event for the engine socket. If it's an event
        // indicating that there's input then tell the user it's okay
        // to read from the engine socket.
        LIBSIM_API_LEAVE1(VisItDetectInput,
                          "WAIT_OBJECT_0+1: Engine socket input. return %d",
                          2);
        return 2;
    }
    else if(n == WAIT_OBJECT_0+2)
    {
        // we received a stdin event. Tell the client that it's
        // okay to read from the console.
        LIBSIM_API_LEAVE1(VisItDetectInput,
                          "WAIT_OBJECT_0+2: Console socket input. return %d",
                          3);
        fprintf(stderr, "WAIT_OBJECT_0+2\n");
        return 3;
    }

    LIBSIM_API_LEAVE1(VisItDetectInput,
                      "Unspecified error. return %d",
                      -5);
    return -5;
}
#else
int
VisItDetectInputWithTimeout(int blocking, int timeoutVal,
    int consoleFileDescriptor)
{
   /*  RETURN CODES:
      -5: Logic error (fell through all cases)
      -4: Logic error (no descriptors but blocking)
      -3: Logic error (a socket was selected but not one we set)
      -2: Unknown error in select
      -1: Interrupted by EINTR in select
       0: Okay - Timed out
       1: Listen  socket input
       2: Engine  socket input
       3: Console socket input
   */
   int maxDescriptor = MAX(MAX(consoleFileDescriptor,
                               engineSocket),
                               listenSocket);

   fd_set readSet;
   int status = 0;
   struct timeval TimeoutValue = {0,0};
   struct timeval *timeout = (blocking ? NULL : &TimeoutValue);
   TimeoutValue.tv_sec = timeoutVal / 1000000;
   TimeoutValue.tv_usec = timeoutVal - (TimeoutValue.tv_sec * 1000000);

   LIBSIM_API_ENTER2(VisItDetectInput, "blocking=%d, consoleFile=%d",
                     blocking, consoleFileDescriptor);

   if (consoleFileDescriptor < 0 &&
       engineSocket <= VISIT_INVALID_SOCKET &&
       listenSocket <= VISIT_INVALID_SOCKET)
   {
      if (blocking)
      {
         LIBSIM_API_LEAVE1(VisItDetectInput,
                           "Logic error (no descriptors but blocking). return %d",
                           -4);
         return -4;
      }
      else
      {
         LIBSIM_API_LEAVE1(VisItDetectInput,
                           "Okay - Timed out. return %d",
                           0);
         return 0;
      }
   }

   FD_ZERO(&readSet);

   /* If we were told to, select on the console socket */
   if (consoleFileDescriptor >= 0)
      FD_SET(consoleFileDescriptor, &readSet);

   /* If we're connected, select on the control socket */
   if (engineSocket >= 0)
      FD_SET(engineSocket, &readSet);

   /* If we're connected, do *not* select on the listen socket */
   /* This forces us to have only one client at a time. */
   if (engineSocket < 0 &&
       listenSocket >= 0)
      FD_SET(listenSocket, &readSet);

   status = select(maxDescriptor+1,
                   &readSet, (fd_set*)NULL, (fd_set*)NULL, timeout);

   if (status < 0)
   {
      if (errno == EINTR)
      {
         /* Interruption will likely cause a program exit anyway */
         LIBSIM_API_LEAVE1(VisItDetectInput,
                           "Interrupted by EINTR in select. return %d",
                           -1);
         return -1;
      }
      else
      {
         /* This should never happen; internal error */
         LIBSIM_API_LEAVE1(VisItDetectInput,
                           "Unknown error in select. return %d",
                           -2);
         return -2;
      }
   }
   else if (status == 0)
   {
      /* We timed out */
      LIBSIM_API_LEAVE1(VisItDetectInput,
                        "Okay - Timed out. return %d",
                        0);
      return 0;
   }
   else
   {
      /* Determine if it's a new connection attempt or
         commands from an existing one */
      if (listenSocket >= 0 &&
          FD_ISSET(listenSocket, &readSet))
      {
         LIBSIM_API_LEAVE1(VisItDetectInput,
                           "Listen  socket input. return %d",
                           1);
         return 1;
      }
      else if (engineSocket >= 0 &&
               FD_ISSET(engineSocket, &readSet))
      {
         LIBSIM_API_LEAVE1(VisItDetectInput,
                           "Engine  socket input. return %d",
                           2);
         return 2;
      }
      else if (consoleFileDescriptor >= 0 &&
               FD_ISSET(consoleFileDescriptor, &readSet))
      {
         LIBSIM_API_LEAVE1(VisItDetectInput,
                           "Console socket input. return %d",
                           3);
         return 3;
      }
      else
      {
         LIBSIM_API_LEAVE1(VisItDetectInput,
                           "Logic error (a socket was selected but not one we set). return %d",
                           -3);
         return -3;
      }
   }
   /*return -5; Compilers complain because this line cannot be reached */
}
#endif

/*******************************************************************************
*
* Name: VisItAttemptToCompleteConnection
*
* Purpose: Accept the socket, verify security keys, get the connection
*          parameters from the client, load the VisIt engine library,
*          create the Engine and connect to the viewer.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Fri Jul 25 15:54:14 PDT 2008
*   Trace information.
*
*   Brad Whitlock, Thu Mar 26 13:38:50 PDT 2009
*   Install our command callback handler.
*
*   Brad Whitlock, Fri Nov 26 00:23:34 PDT 2010
*   Add Windows code
*
*******************************************************************************/
int VisItAttemptToCompleteConnection(void)
{
    VISIT_SOCKET socket;

    LIBSIM_API_ENTER(VisItAttemptToCompleteConnection);

    /* wait for a connection -- only process 0 does this */
    if (parallelRank == 0)
    {
        socket = AcceptConnection();

        if (socket < 0)
        {
            LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, 
                              "socket<0, return %d", VISIT_ERROR);
            return VISIT_ERROR;
        }
    }

    /* verify security keys */
    if (!VerifySecurityKeys(socket))
    {
        LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, 
                          "VerifySecurityKeys failed. return %d", VISIT_ERROR);
        return VISIT_ERROR;
    }

    /* get the connection parameters */
    if (!GetConnectionParameters(socket))
    {
        LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, 
                          "GetConnectionParameters failed. return %d", VISIT_ERROR);
        return VISIT_ERROR;
    }

    /* load the library */
    if (LoadVisItLibrary() == 0)
    {
        LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, 
                          "LoadVisItLibrary failed. return %d", VISIT_ERROR);
        return VISIT_ERROR;
    }

    /* connect to the viewer */
    if (CreateEngineAndConnectToViewer() == 0)
    {
        LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, 
                          "CreateEngineAndConnectToViewer failed. return %d", VISIT_ERROR);
        return VISIT_ERROR;
    }

    /* get the socket for listening from the viewer */
    if (parallelRank == 0)
    {
        LIBSIM_MESSAGE("Calling visit_getdescriptor");
        engineSocket = callbacks->control.get_descriptor(engine);
        LIBSIM_MESSAGE1("visit_getdescriptor returned %d", (int)engineSocket);
#ifdef _WIN32
#ifndef SELECT_ENGINE_SOCKET
        WSAEventSelect(engineSocket, engineevent, FD_READ);
#endif
        SetEvent(listeneventCB);
#endif
    }

    /* Install our local command callback handler so we can monitor the
     * callbacks before handing them off to the user's command callback. This
     * lets us implement synchronization for the user.
     */
    if(callbacks->control.set_command_callback != NULL)
        (*callbacks->control.set_command_callback)(engine, visit_handle_command_callback, NULL);

    LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, "return %d", VISIT_OKAY);
    return VISIT_OKAY;
}

/*******************************************************************************
*
* Name: VisItSetSlaveProcessCallback
*
* Purpose: Set the callback to inform slave processes that they should
*          call VisItProcessEngineCommand.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*   Brad Whitlock, Fri Jul 25 15:55:53 PDT 2008
*   Trace information.
*
*******************************************************************************/
void VisItSetSlaveProcessCallback(void (*spic)(void))
{
    LIBSIM_API_ENTER1(VisItSetSlaveProcessCallback, "spic=%p", (void*)spic);
    LIBSIM_MESSAGE("Calling visit_set_slave_process_callback");
    if(callbacks != NULL && callbacks->control.set_slave_process_callback)
        (*callbacks->control.set_slave_process_callback)(spic);
    LIBSIM_API_LEAVE(VisItSetSlaveProcessCallback);
}

/*******************************************************************************
*
* Name: VisItSetCommandCallback
*
* Purpose: Set the callback for processing control commands.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 15:56:41 PDT 2008
*  Trace information.
*
*******************************************************************************/
void VisItSetCommandCallback(void (*scc)(const char*,const char*,void*),
    void *sccdata)
{
    LIBSIM_API_ENTER1(VisItSetCommandCallback, "scc=%p", (void*)scc);
    LIBSIM_MESSAGE("Calling visit_set_command_callback");

    visit_command_callback = scc;
    visit_command_callback_data = sccdata;

    LIBSIM_API_LEAVE(VisItSetCommandCallback);
}

/*******************************************************************************
*
* Name: VisItProcessEngineCommand
*
* Purpose: Process requests coming from the client.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 15:56:41 PDT 2008
*  Trace information.
*
*  Brad Whitlock, Fri Nov 26 00:25:24 PDT 2010
*  Add Windows code.
*
*******************************************************************************/
int VisItProcessEngineCommand(void)
{
    int retval = VISIT_ERROR;
    LIBSIM_API_ENTER(VisItProcessEngineCommand);
    if (!engine)
    {
        LIBSIM_MESSAGE("engine == NULL");
        retval = VISIT_OKAY;
    }
    else if (callbacks != NULL)
    {
        LIBSIM_MESSAGE("Calling visit_processinput");
        retval = ((*callbacks->control.process_input)(engine) == 1) ?
            VISIT_OKAY : VISIT_ERROR;
        LIBSIM_MESSAGE1("visit_processinput returned: %d", retval);
    }
#ifdef _WIN32
    if(parallelRank == 0)
        SetEvent(engineeventCB);
#endif
    LIBSIM_API_LEAVE1(VisItProcessEngineCommand, "return %d", retval);
    return retval;
}

/*******************************************************************************
*
* Name: VisItTimeStepChanged
*
* Purpose: Tell VisIt a new time step is ready.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 16:00:32 PDT 2008
*  Trace information.
*
*******************************************************************************/
void VisItTimeStepChanged(void)
{
    LIBSIM_API_ENTER(VisItTimeStepChanged);
    /* Make sure the function exists before using it. */
    if (engine && callbacks != NULL && callbacks->control.time_step_changed)
    {
        LIBSIM_MESSAGE("Calling visit_time_step_changed");
        (*callbacks->control.time_step_changed)(engine);
    }
    LIBSIM_API_LEAVE(VisItTimeStepChanged);
}

/*******************************************************************************
*
* Name: VisItUpdatePlots
*
* Purpose: Tell VisIt to update its plots.
*
* Author: Brad Whitlock, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 16:00:32 PDT 2008
*  Trace information.
*
*******************************************************************************/
void VisItUpdatePlots(void)
{
    LIBSIM_API_ENTER(VisItUpdatePlots);

    /* Make sure the function exists before using it. */
    if (engine && callbacks != NULL && callbacks->control.execute_command)
    {
        LIBSIM_MESSAGE("Calling visit_execute_command: UpdatePlots");
        (*callbacks->control.execute_command)(engine, "UpdatePlots");

        if(visit_sync_enabled)
            VisItSynchronize();
    }
    LIBSIM_API_LEAVE(VisItUpdatePlots);
}

/*******************************************************************************
*
* Name: VisItExecuteCommand
*
* Purpose: Tell VisIt to execute a command.
*
* Author: Brad Whitlock, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 16:00:32 PDT 2008
*  Trace information.
*
*******************************************************************************/
void VisItExecuteCommand(const char *command)
{
    LIBSIM_API_ENTER(VisItExecuteCommand);
    /* Make sure the function exists before using it. */
    if (engine && callbacks != NULL && callbacks->control.execute_command)
    {
        char *cmd2 = NULL;
        LIBSIM_MESSAGE("Calling visit_execute_command");
        cmd2 = (char *)malloc(strlen(command) + 1 + 10);
        sprintf(cmd2, "Interpret:%s", command);
        (*callbacks->control.execute_command)(engine, cmd2);
        free(cmd2);

        if(visit_sync_enabled)
            VisItSynchronize();
    }
    LIBSIM_API_LEAVE(VisItExecuteCommand);
}

/*******************************************************************************
*
* Name: VisItDisconnect
*
* Purpose: Disconnect from the client.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 16:00:32 PDT 2008
*  Trace information.
*
*  Brad Whitlock, Fri Feb 13 09:45:14 PST 2009
*  Delete the callbacks structure and set it to NULL.
*
*******************************************************************************/
void VisItDisconnect(void)
{
    LIBSIM_API_ENTER(VisItDisconnect);
    LIBSIM_MESSAGE("Calling visit_disconnect");
    if(callbacks != NULL)
    {
        if(callbacks->control.disconnect)
            (*callbacks->control.disconnect)();

        free(callbacks);
        callbacks = NULL;
    }
    engineSocket = VISIT_INVALID_SOCKET;
    engine = NULL;
#ifdef _WIN32
    inputThreadsStarted = 0;
#endif
    CloseVisItLibrary();
    LIBSIM_API_LEAVE(VisItDisconnect);
}

/*******************************************************************************
*
* Name: VisItIsConnected
*
* Purpose: Returns whether VisIt is connected
*
* Author: Brad Whitlock, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*
*******************************************************************************/
int VisItIsConnected(void)
{
    return ((engine != 0) && (callbacks != NULL)) ? 1 : 0;
}

/*******************************************************************************
*
* Name: VisItSaveImage
*
* Purpose: Tell VisIt to save the active network as an image.
*
* Author: Brad Whitlock, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*
*******************************************************************************/
int VisItSaveWindow(const char *filename, int w, int h, int format)
{
    int ret = VISIT_ERROR;
    LIBSIM_API_ENTER(VisItSaveWindow);
    /* Make sure the function exists before using it. */
    if (engine && callbacks != NULL && callbacks->control.save_window)
    {
        LIBSIM_MESSAGE("Calling visit_save_window");
        ret = (*callbacks->control.save_window)(engine, filename, w, h, format);
    }
    LIBSIM_API_LEAVE(VisItSaveWindow);
    return ret;
}

/*******************************************************************************
*
* Name: VisItGetLastError
*
* Purpose: Returns the last error message generated.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 16:00:32 PDT 2008
*  Trace information.
*
*******************************************************************************/
char *VisItGetLastError()
{
    char *err = NULL;

    LIBSIM_API_ENTER(VisItGetLastError);
    err = strdup(lastError);
    strcpy(lastError, "");
    LIBSIM_API_LEAVE1(VisItGetLastError, "return %s", err);
    return err;
}

/******************************************************************************
*
* Name: VisItOpenTraceFile
*
* Purpose: Opens the trace file that libsim will use.
*
* Programmer: Brad Whitlock
* Date:       Mon Feb  9 10:01:08 PST 2009
*
* Modifications:
*
******************************************************************************/
void VisItOpenTraceFile(const char *filename)
{
    simv2_set_trace_file(fopen(filename, "wt"));
}

/******************************************************************************
*
* Name: VisItCloseTraceFile
*
* Purpose: Opens the trace file that libsim will use.
*
* Programmer: Brad Whitlock
* Date:       Mon Feb  9 10:01:08 PST 2009
*
* Modifications:
*
******************************************************************************/
void VisItCloseTraceFile(void)
{
    simv2_set_trace_file(NULL);
}

/******************************************************************************
*
* Name: VisItSet*
*
* Purpose: Define the bodies of all of the functions that let the user
*          register their libsim callback function.
*
* Programmer: Brad Whitlock
* Date:       Mon Feb  9 10:01:08 PST 2009
*
* Modifications:
*
******************************************************************************/

#define VISIT_SET_CALLBACK_BODY(F, T) \
int VisItSet##F(T, void *cbdata) \
{ \
    int retval = VISIT_ERROR; \
    LIBSIM_API_ENTER(VisIt##F);\
    if (engine && callbacks != NULL && callbacks->data.set_##F)\
    {\
        LIBSIM_MESSAGE("Calling VisIt"#F);\
        (*callbacks->data.set_##F)(cb, cbdata);\
        retval = VISIT_OKAY; \
    }\
    else\
        fprintf(stderr, "VisIt"#F" should not be called until VisIt connects to the simulation.\n");\
    LIBSIM_API_LEAVE(VisIt##F); \
    return retval;\
}

DECLARE_DATA_CALLBACKS(VISIT_SET_CALLBACK_BODY, CB_ARGS)

/******************************************************************************
*
* Name: VisItSynchronize
*
* Purpose: Sends a synchronize message to VisIt if we're connected to VisIt
*          and spawns a sub-eventloop to wait for the sync message to be
*          returned from VisIt. This permits us to write code in the sim that
*          blocks until it's done and still processes events from VisIt needed
*          to get us to the point where the message gets read.
*
* Note:    This function has not been tried for sims that use polling.
*
* Programmer: Brad Whitlock
* Date:       Thu Mar 26 13:28:11 PDT 2009
*
* Modifications:
*
******************************************************************************/

static void
visit_sync_helper(void *cbdata)
{
    int *syncing = (int *)cbdata;
    *syncing = 0;
}

int
VisItSynchronize(void)
{
    int blocking = 1;
    int syncing = 1;
    int visitstate = 0, err = 0;

    LIBSIM_API_ENTER(VisItSynchronize);

    if(!VisItIsConnected())
    {
        return VISIT_OKAY;
    }

    /* Send a sync to the viewer. When we get it back the loop will end. */
    visit_add_sync(visit_sync_helper, &syncing);

    do
    {
        /* Get input from VisIt. */
        if(parallelRank == 0)
            visitstate = VisItDetectInput(blocking, -1);
        BroadcastInt(&visitstate, 0);

        /* Do different things depending on the output from VisItDetectInput. */
        if(visitstate >= -5 && visitstate <= -1)
        {
            fprintf(stderr, "Can't recover from error!\n");
            err = 1;
        }
        else if(visitstate == 0)
        {
            /* There was no input from VisIt. We're blocking so this should not happen. */
        }
        else if(visitstate == 1)
        {
            /* VisIt is trying to connect to the sim. We're already connected 
             * so this can't happen.
             */
        }
        else if(visitstate == 2)
        {
            /* VisIt wants to tell the engine something. */
            if(!visit_process_engine_command())
            {
                /* Disconnect on an error or closed connection. */
                VisItDisconnect();
                err = 1;
            }
        }
        else if(visitstate == 3)
        {
            /* We're not trapping for console input. */
        }
    } while(syncing && err == 0);

    LIBSIM_API_LEAVE(VisItSynchronize);
    return (err==0) ? VISIT_OKAY : VISIT_ERROR;
}

/******************************************************************************
*
* Name: VisItEnableSynchronize
*
* Purpose: This method lets us turn off automatic synchronization for functions
*          that communicate with VisIt's viewer.
*
* Programmer: Brad Whitlock
* Date:       Thu Mar 26 13:28:11 PDT 2009
*
* Modifications:
*
******************************************************************************/

void
VisItEnableSynchronize(int mode)
{
    visit_sync_enabled = mode ? 1 : 0;
}

/******************************************************************************
*
* Name: VisItDebug*
*
* Purpose: Let the user write to the VisIt debug logs.
*
* Programmer: Brad Whitlock
* Date:       Thu Mar 26 13:28:11 PDT 2009
*
* Modifications:
*
******************************************************************************/

#define DBG_BUFFER_SIZE 1000
#define VISIT_DEBUG_FUNCTION(N)\
void \
VisItDebug##N(const char *format, ...)\
{\
    char buffer[DBG_BUFFER_SIZE]; \
    va_list args; \
    va_start(args, format); \
    vsnprintf(buffer, DBG_BUFFER_SIZE, format, args);\
    va_end(args);\
    if(engine && callbacks != NULL && callbacks->control.debug_logs != NULL)\
        (*callbacks->control.debug_logs)(N, buffer);\
}

VISIT_DEBUG_FUNCTION(1)
VISIT_DEBUG_FUNCTION(2)
VISIT_DEBUG_FUNCTION(3)
VISIT_DEBUG_FUNCTION(4)
VISIT_DEBUG_FUNCTION(5)
