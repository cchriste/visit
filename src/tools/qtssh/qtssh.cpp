#include <ViewerPasswordWindow.h>
#include <qapplication.h>
#include <windows.h>
#include <stdlib.h>
#include <io.h>

// Include the header file for the remote command SSH library.
#include <RemoteCommand.h>

// Globals
static char *username = 0;

// ****************************************************************************
// Function: graphicalGetPassword
//
// Purpose: This is the callback function for getting the password from the
//          user. It uses the ViewerPasswordWindow class to get the password
//          graphically.
//
// Programmer: Brad Whitlock
// Creation:   Thu Aug 29 14:57:24 PST 2002
//
// Modifications:
//   Brad Whitlock, Mon Feb 23 15:03:24 PST 2004
//   I made the okay parameter be 0 if the password window returns NULL
//   for a password. This lets us quit without crashing.
//
// ****************************************************************************

const char *
graphicalGetPassword(const char *host, int *okay)
{
    const char *retval = ViewerPasswordWindow::getPassword(username, host);
    *okay = (retval != 0);
    return retval;
}

// ****************************************************************************
// Function: main
//
// Purpose:
//   This is the main function for the qtssh program.
//
// Programmer: Brad Whitlock
// Creation:   Thu Aug 29 14:59:21 PST 2002
//
// Modifications:
//   Brad Whitlock, Fri Oct 10 14:24:27 PST 2003
//   Added the -p argument.
//
//   Brad Whitlock, Tue Dec 21 16:17:12 PST 2004
//   Added code to close stdin so we don't gobble up input typed into the CLI.
//
//   Jeremy Meredith, Wed Jun 27 12:08:52 EDT 2007
//   Added support for standard SSH-style remote port forwarding arguments
//   of the form "-R rp:lh:lp:".
//
// ****************************************************************************

int
main(int argc, char *argv[])
{
    const char *host = "localhost";
    bool userSpecified = false;
    bool first = true;
    const char **commands = new const char*[100];
    int i, command_count = 0;
    bool printArgs = false;
    bool hostSpecified = false;
    int  port = 22;
    char portfwd[1024] = "";
    int portfwdpos = 0;

    // Initialize the command line array.
    for(i = 0; i < 100; ++i)
        commands[i] = 0;

    // Create a new application.
    new QApplication(argc, argv);

    // Get the username.
    DWORD namelen = 100;
    username = new char[100];
    GetUserName(username, &namelen);
  
    // Parse the command line arguments that matter.
    for(i = 1; i < argc; ++i)
    {
        char *arg = argv[i];
        if(arg[0] == '-')
        {
            if(strcmp(arg, "-l") == 0)
            {
                if(i+1 < argc)
                {
                    delete [] username;
                    username = argv[i+1];
                    ++i;
                }
            }
            else if(strcmp(arg, "-p") == 0)
            {
                if(i+1 < argc)
                {
                    int tempPort = atoi(argv[i+1]);
                    if(tempPort >= 0)
                        port = tempPort;
                    ++i;
                }
            }
            else if(strcmp(arg, "-R") == 0)
            {
                // we must compose the port forward string
                // according to PUTTY.H
                if(i+1 < argc)
                {
                    int   l = strlen(argv[i+1]);
                    char *s = strdup(argv[i+1]);
                    char *p = s;
                    int part = 0;
                    char *remoteport;
                    char *localhost;
                    char *localport;
                    for (int j=0; j<l+1; j++)
                    {
                        if (s[j] == ':' || s[j] == '\0')
                        {
                            s[j] = '\0';
                            switch (part)
                            {
                              case 0:
                                remoteport = p;
                                break;
                              case 1:
                                localhost = p;
                                break;
                              case 2:
                                localport = p;
                                break;
                              default:
                                // we'll detect this as an error automatically
                                break;
                            }
                            p = &(s[j+1]);
                            part++;
                        }
                    }
                    if (part == 3)
                    {
                        // first is the R, for Remote
                        portfwd[portfwdpos++] = 'R';
                        for (int j=0; j<strlen(remoteport); j++)
                            portfwd[portfwdpos++] = remoteport[j];
                        // remote/local separated by tab
                        portfwd[portfwdpos++] = '\t';
                        for (int j=0; j<strlen(localhost); j++)
                            portfwd[portfwdpos++] = localhost[j];
                        // host/port separated by a colon
                        portfwd[portfwdpos++] = ':';
                        for (int j=0; j<strlen(localport); j++)
                            portfwd[portfwdpos++] = localport[j];
                        // terminate this forward with a \0
                        portfwd[portfwdpos++] = '\0';
                        // and the final forward is terminated with another \0
                        portfwd[portfwdpos] = '\0';
                    }
                    ++i;
                }
            }
            else if(strcmp(arg, "-D") == 0)
            {
                printArgs = true;
            }
            else if(!first)
            {
                // Store the pointer into the option list.
                commands[command_count++] = arg;
            }
        }
        else if(first)
        {
            // The first option without a '-' character must be the hostname.
            first = false;
            hostSpecified = true;
            host = arg;
        }
        else
        {
            // Store the pointer into the option list.
            commands[command_count++] = arg;
        }
    }

    if(printArgs)
    {
        qDebug("username=%s", username);
        qDebug("host=%s", host);
        for(int j = 0; j < command_count; ++j)
            qDebug("command[%d]=%s", j, commands[j]);
    }

    if(!hostSpecified)
    {
        qDebug("usage: %s [options] hostname [command args...]", argv[0]);
        qDebug("");
        qDebug("options:");
        qDebug("    -l username     Sets the user login name.");
        qDebug("    -p portnum      Sets the port number used to connect.");
        qDebug("    -D              Prints command line arguments.");
        return -1;
    }

    // Close stdin so we don't try and intercept input from the CLI.
    _close(0);

    // Run the command on the remote machine.
    RunRemoteCommand(username, host, port, commands, command_count,
                     graphicalGetPassword, 1, portfwd);

    // Clean up.
    if(!userSpecified)
        delete [] username;
    delete [] commands;

    return 0;
}
