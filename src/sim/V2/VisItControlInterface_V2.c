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

#include "VisItControlInterface_V1.h"
#include "SimV2Tracing.h"

#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <netdb.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <pwd.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#include <visit-config.h> /* For HAVE_SOCKLEN_T */

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


/* VisIt Engine Library function pointers */
typedef struct
{
    void *(*v_getengine)(void);
    int   (*v_getdescriptor)(void*);
    int   (*v_processinput)(void*);
    int   (*v_initialize)(void*,int,char**);
    int   (*v_connectviewer)(void*,int,char**);
    void  (*v_time_step_changed)(void*);
    void  (*v_update_plots)(void*);
    void  (*v_execute_command)(void*,const char*);
    void  (*v_disconnect)();
    void  (*v_set_slave_process_callback)(void(*)());
    void  (*v_set_command_callback)(void*,void(*)(const char*,int,
                                                  float,const char*));
} control_callback_t;

static control_callback_t *callbacks = NULL;


/* Internal Variables */
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
static int         listenSocket = -1;
struct sockaddr_in listenSockAddr;
static int         engineSocket = -1;
static int         isParallel = FALSE;
static int         parallelRank = 0;
static void       *dl_handle;
static char        lastError[1024] = "";


/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
                               Internal Functions
 *******************************************************************************
 *******************************************************************************
 ******************************************************************************/



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
static void ReceiveSingleLineFromSocket(char *buffer, size_t maxlen, int desc)
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
static char *ReceiveContinuousLineFromSocket(char *buffer, size_t maxlen, int desc)
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
static int SendStringOverSocket(char *buffer, int desc)
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
static int VerifySecurityKeys(int desc)
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
static int GetConnectionParameters(int desc)
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
    LIBSIM_MESSAGE("Calling v_engine");
    engine = (*callbacks->v_getengine)();
    if (!engine)
    {
        LIBSIM_API_LEAVE1(CreateEngineAndConnectToViewer,
                         "engine could not be allocated. return %d",
                         FALSE);
        return FALSE;
    }

    LIBSIM_MESSAGE_STRINGLIST("Calling v_initialize: argv",
                              engine_argc, engine_argv);
    if (!(*callbacks->v_initialize)(engine, engine_argc, engine_argv))
    {
        VisItDisconnect();
        LIBSIM_API_LEAVE1(CreateEngineAndConnectToViewer,
                         "v_initialize failed. return %d",
                         FALSE);
        return FALSE;
    }

    LIBSIM_MESSAGE_STRINGLIST("Calling v_connectviewer: argv",
                              engine_argc, engine_argv);
    if (!(*callbacks->v_connectviewer)(engine, engine_argc, engine_argv))
    {
        VisItDisconnect();
        LIBSIM_API_LEAVE1(CreateEngineAndConnectToViewer,
                         "v_connectviewer failed. return %d", 
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

    LIBSIM_MESSAGE("Calling gethostbyname");
    localhostEnt = gethostbyname(localhostStr);
    if (localhostEnt == NULL)
    {
        /* Couldn't get the full host entry; it's probably invalid */
        LIBSIM_API_LEAVE1(GetLocalhostName, 
                         "gethostbyname failed. return %d", FALSE);
        return FALSE;
    }
    sprintf(localhost, localhostEnt->h_name);
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
       listenSocket = -1;
       LIBSIM_API_LEAVE1(StartListening, "port not found. return %d", FALSE);
       return FALSE;
    }

    LIBSIM_MESSAGE("Calling listen() on socket");
    err = listen(listenSocket, 5);
    if (err)
    {
       listenSocket = -1;
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
static int AcceptConnection(void)
{
    int desc = -1;
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
    while (desc == -1 && errno != EWOULDBLOCK);

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
    struct passwd *users_passwd_entry = NULL;

    LIBSIM_API_ENTER(GetHomeDirectory);
    users_passwd_entry = getpwuid(getuid());
    LIBSIM_API_LEAVE1(GetHomeDirectory, "homedir=%s", 
                     users_passwd_entry->pw_dir);

    return users_passwd_entry->pw_dir;
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

    snprintf(str, 1024, "%s/.visit", GetHomeDirectory());
    mkdir(str, 7*64 + 7*8 + 7);
    LIBSIM_MESSAGE1("mkdir %s", str);

    snprintf(str, 1024, "%s/.visit/simulations", GetHomeDirectory());
    mkdir(str, 7*64 + 7*8 + 7);
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
* Name: dlsym_function
*
* Purpose: Call dlsym and cast the result to a function pointer.
*          Basically, this exists because the OSF compiler won't let you
*          cast from void* to void(*)() without a warning, and I can at
*          least contain it to a single warning by creating this function.
*
* Author: Jeremy Meredith, B Division, Lawrence Livermore National Laboratory
*
* Modifications:
*  Brad Whitlock, Fri Jul 25 15:00:20 PDT 2008
*  Trace information.
*
*******************************************************************************/

typedef void (*function_pointer)();
static function_pointer dlsym_function(void *h, const char *n)
{
    function_pointer f = NULL;

    LIBSIM_API_ENTER2(dlsym_function, "handle=%p, name=%s", h, n);
    f = (function_pointer)dlsym(h,n);
    LIBSIM_API_LEAVE1(dlsym_function, "func=%p", (void*)f);
    return f;
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
*******************************************************************************/
static int LoadVisItLibrary(void)
{
   char lib[256];

   /* For a while, something was hanging if you tried to re-use the same
      engine library for a second connection.  Things seem better now, so
      there's no need to enable this next line of code.  It is left in as
      a reminder that we may run into the same issue on other platforms.

    if (dl_handle)
        return;
   */

   /* load library */
#ifdef __APPLE__
   const char *extension = "dylib";
#else
   const char *extension = "so";
#endif
   LIBSIM_API_ENTER(LoadVisItLibrary);

   callbacks = (control_callback_t *)malloc(sizeof(control_callback_t));
   memset(callbacks, 0, sizeof(control_callback_t));

   if (isParallel)
   {
       sprintf(lib, "libvisitenginev1_par.%s", extension);
   }
   else
   {
       sprintf(lib, "libvisitenginev1_ser.%s", extension);
   }

   LIBSIM_MESSAGE1("Calling dlopen(%s)", lib);
   dl_handle = dlopen(lib, RTLD_NOW | RTLD_GLOBAL);

   if (!dl_handle && getenv("LD_LIBRARY_PATH"))
   {
      char *libpath = strdup(getenv("LD_LIBRARY_PATH"));
      char *ptr = strstr(libpath, ":");
      if (ptr)
      {
         *ptr = 0;

         if (isParallel) 
             sprintf(lib, "%s/libvisitenginev1_par.%s", libpath, extension);
         else
             sprintf(lib, "%s/libvisitenginev1_ser.%s", libpath, extension);

         LIBSIM_MESSAGE1("Calling dlopen(%s)", lib);
         dl_handle = dlopen(lib, RTLD_NOW | RTLD_GLOBAL);
      }
   }

   if (!dl_handle)
   {
      sprintf(lastError, "Failed to open the VisIt library: %s\n", dlerror());
      LIBSIM_API_LEAVE2(LoadVisItLibrary, "%s: return %d", lastError, FALSE);
      return FALSE;
   }

#define SAFE_DLSYM(f,t,n)                \
   callbacks->f = (t)dlsym_function(dl_handle, n); \
   if (!callbacks->f) \
   { \
      sprintf(lastError, "Failed to open the VisIt library: "\
              "couldn't find symbol '%s': %s\n", n, dlerror()); \
      dl_handle = NULL; LIBSIM_API_LEAVE2(LoadVisItLibrary, "%s: return %d", lastError, FALSE); \
      return FALSE; \
   }

   SAFE_DLSYM(v_getengine, void *(*)(void), "get_engine");
   SAFE_DLSYM(v_getdescriptor, int (*)(void *), "get_descriptor");
   SAFE_DLSYM(v_processinput, int (*)(void *), "process_input");
   SAFE_DLSYM(v_initialize, int (*)(void *, int, char **), "initialize");
   SAFE_DLSYM(v_connectviewer, int (*)(void *, int, char **), "connect_to_viewer");
   SAFE_DLSYM(v_time_step_changed, void (*)(void *), "time_step_changed");
   SAFE_DLSYM(v_update_plots, void (*)(void *), "update_plots");
   SAFE_DLSYM(v_execute_command, void (*)(void *,const char*), "execute_command");
   SAFE_DLSYM(v_disconnect, void (*)(), "disconnect");
   SAFE_DLSYM(v_set_slave_process_callback, void (*)(void (*)()), "set_slave_process_callback");
   SAFE_DLSYM(v_set_command_callback, void (*)(void*,void (*)(const char*,int,float,const char*)), "set_command_callback");

#ifdef VISIT_QUERIES_VIA_LIBSIM
   /* Load the VisIt query functions.*/
   LoadVisItQueries();
#endif

   LIBSIM_API_LEAVE1(LoadVisItLibrary, "return %d", TRUE);
   return TRUE;
}

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
   snprintf(command, 1024, "%s -compiler %s %s -env -engine 2>/dev/null",
           visitpath, STR(VISIT_COMPILER), visit_options ? visit_options : "");
#else
   snprintf(command, 1024, "%s %s -env -engine 2>/dev/null",
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

   LIBSIM_API_LEAVE1(ReadEnvironmentFromCommand, "return %d", (ptr-output));
   return (ptr - output);
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
*******************************************************************************/
int VisItSetupEnvironment(void)
{
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
    FILE *file;
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
        snprintf(simulationFileName, 255, "%s/.visit/simulations/%012d.%s.sim2",
                 GetHomeDirectory(), (int)time(NULL), name);
    }
    else
    {
         snprintf(simulationFileName, 255, "%s", absoluteFilename);
    }

    LIBSIM_MESSAGE1("Opening sime file %s", simulationFileName);  
    file = fopen(simulationFileName, "wt");
    if (!file)
    {
        LIBSIM_API_LEAVE1(VisItInitializeSocketAndDumpSimFile, "return %d", FALSE);
        return FALSE;
    }

    atexit(RemoveSimFile);

    if (!GetLocalhostName())
    {
        fclose(file);
        LIBSIM_API_LEAVE1(VisItInitializeSocketAndDumpSimFile, "return %d", FALSE);
        return FALSE;
    }

    if (!StartListening())
    {
        fclose(file);
        LIBSIM_API_LEAVE1(VisItInitializeSocketAndDumpSimFile, "return %d", FALSE);
        return FALSE;
    }

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
int  VisItDetectInput(int blocking, int consoleFileDescriptor)
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
   struct timeval ZeroTimeout = {0,0};
   struct timeval *timeout = (blocking ? NULL : &ZeroTimeout);

   LIBSIM_API_ENTER2(VisItDetectInput, "blocking=%d, consoleFile=%d",
                     blocking, consoleFileDescriptor);

   if (consoleFileDescriptor < 0 &&
       engineSocket < 0 &&
       listenSocket < 0)
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
*******************************************************************************/
int VisItAttemptToCompleteConnection(void)
{
    int socket;

    LIBSIM_API_ENTER(VisItAttemptToCompleteConnection);

    /* wait for a connection -- only process 0 does this */
    if (parallelRank == 0)
    {
        socket = AcceptConnection();

        if (socket < 0)
        {
            LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, 
                              "socket<0, return %d", FALSE);
            return FALSE;
        }
    }

    /* verify security keys */
    if (!VerifySecurityKeys(socket))
    {
        LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, 
                          "VerifySecurityKeys failed. return %d", FALSE);
        return FALSE;
    }

    /* get the connection parameters */
    if (!GetConnectionParameters(socket))
    {
        LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, 
                          "GetConnectionParameters failed. return %d", FALSE);
        return FALSE;
    }

    /* load the library */
    if (LoadVisItLibrary() == 0)
    {
        LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, 
                          "LoadVisItLibrary failed. return %d", FALSE);
        return FALSE;
    }

    /* connect to the viewer */
    if (CreateEngineAndConnectToViewer() == 0)
    {
        LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, 
                          "CreateEngineAndConnectToViewer failed. return %d", FALSE);
        return FALSE;
    }

    /* get the socket for listening from the viewer */
    if (parallelRank == 0)
    {
        LIBSIM_MESSAGE("Calling v_getdescriptor");
        engineSocket = callbacks->v_getdescriptor(engine);
        LIBSIM_MESSAGE1("v_getdescriptor returned %d", (int)engineSocket);
    }

    LIBSIM_API_LEAVE1(VisItAttemptToCompleteConnection, "return %d", TRUE);
    return TRUE;
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
void VisItSetSlaveProcessCallback(void (*spic)())
{
    LIBSIM_API_ENTER1(VisItSetSlaveProcessCallback, "spic=%p", (void*)spic);
    LIBSIM_MESSAGE("Calling v_set_slave_process_callback");
    (*callbacks->v_set_slave_process_callback)(spic);
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
void VisItSetCommandCallback(void (*scc)(const char*,int,float,const char*))
{
    LIBSIM_API_ENTER1(VisItSetCommandCallback, "scc=%p", (void*)scc);
    LIBSIM_MESSAGE("Calling v_set_command_callback");
    (*callbacks->v_set_command_callback)(engine,scc);
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
*******************************************************************************/
int VisItProcessEngineCommand(void)
{
    int retval;
    LIBSIM_API_ENTER(VisItProcessEngineCommand);
    if (!engine)
    {
        LIBSIM_MESSAGE("engine == NULL");
        retval = FALSE;
    }
    else
    {
        LIBSIM_MESSAGE("Calling v_processinput");
        retval = (*callbacks->v_processinput)(engine);
    }
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
    if (engine && callbacks->v_time_step_changed)
    {
        LIBSIM_MESSAGE("Calling v_time_step_changed");
        (*callbacks->v_time_step_changed)(engine);
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
    if (engine && callbacks->v_update_plots)
    {
        LIBSIM_MESSAGE("Calling v_update_plots");
        (*callbacks->v_update_plots)(engine);
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
    if (engine && callbacks->v_execute_command)
    {
        LIBSIM_MESSAGE("Calling v_execute_command");
        (*callbacks->v_execute_command)(engine, command);
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
*******************************************************************************/
void VisItDisconnect(void)
{
    LIBSIM_API_ENTER(VisItDisconnect);
    LIBSIM_MESSAGE("Calling v_disconnect");
    (*callbacks->v_disconnect)();
    engineSocket = -1;
    engine = 0;
    LIBSIM_API_LEAVE(VisItDisconnect);
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
 *****************************************************************************/
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
 *****************************************************************************/
void VisItCloseTraceFile(void)
{
    simv2_set_trace_file(NULL);
}
