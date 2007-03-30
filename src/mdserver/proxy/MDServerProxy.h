// ************************************************************************* //
//                              MDServerProxy.h                              //
// ************************************************************************* //

#ifndef _MSSERVERPROXY_H_
#define _MDSERVERPROXY_H_
#include <mdsproxy_exports.h>

#include <RemoteProxyBase.h>
#include <GetDirectoryRPC.h>
#include <ChangeDirectoryRPC.h>
#include <CloseDatabaseRPC.h>
#include <ConnectCallback.h>
#include <ConnectRPC.h>
#include <CreateGroupListRPC.h>
#include <ExpandPathRPC.h>
#include <GetFileListRPC.h>
#include <GetMetaDataRPC.h>
#include <GetSILRPC.h>
#include <avtDatabaseMetaData.h>
#include <SILAttributes.h>

#include <string>
#include <vector>
#include <maptypes.h>
#include <vectortypes.h>

// ****************************************************************************
//  Class: MDServerProxy
//
//  Purpose:
//      MDServerProxy is a proxy class for creating and controlling a viewer.
//
//  Note:
//
//  Programmer: Jeremy Meredith
//  Creation:   August 18, 2000
//
//  Modifications:
//    Brad Whitlock, Mon Aug 21 14:51:36 PST 2000
//    I added a cc, assignment operator to FileList.
//
//    Jeremy Meredith, Fri Sep  1 17:17:12 PDT 2000
//    Added GetMetaData() and its associated objects (RPC, data).
//
//    Eric Brugger, Wed Oct 25 15:30:02 PDT 2000
//    I removed the argument "prog" from the Create method.
//
//    Brad Whitlock, Thu Nov 16 16:21:43 PST 2000
//    I added code that allows it to connect to a process that did not
//    launch it.
//
//    Sean Ahern, Wed Feb 28 14:34:26 PST 2001
//    Added the CreateGroupListRPC.
//
//    Hank Childs, Thu Mar 29 16:49:59 PST 2001
//    Added GetSIL.
//
//    Brad Whitlock, Tue Feb 12 13:56:43 PST 2002
//    Added ExpandPath.
//
//    Brad Whitlock, Mon Mar 25 16:17:24 PST 2002
//    Took away some methods.
//
//    Brad Whitlock, Fri May 3 17:22:54 PST 2002
//    Added a method for setting the remote user name.
//
//    Brad Whitlock, Tue Jul 30 10:52:33 PDT 2002
//    I added a method to close the active database.
//
//    Brad Whitlock, Mon Aug 26 15:53:27 PST 2002
//    I added methods to get the file separator.
//
//    Brad Whitlock, Fri Sep 27 15:34:11 PST 2002
//    I added methods to support a launch progress window.
//
//    Brad Whitlock, Thu Dec 26 16:17:27 PST 2002
//    I added support for security keys.
//
//    Brad Whitlock, Mon Mar 24 11:55:35 PDT 2003
//    I made VisIt capable of automatically detecting time-varying databases.
//
//    Brad Whitlock, Tue Apr 22 13:50:40 PST 2003
//    I fixed a crash on Windows.
//
//    Brad Whitlock, Fri May 2 12:14:28 PDT 2003
//    I made it inherit from RemoteProxyBase.
//
//    Brad Whitlock, Tue May 13 15:29:50 PST 2003
//    I made it possible to get metadata for a certain timestate.
//
// ****************************************************************************

class MDSERVER_PROXY_API MDServerProxy : public RemoteProxyBase
{
public:
    struct FileEntry
    {
        FileEntry();
        FileEntry(const FileEntry &);
        ~FileEntry();
        void operator = (const FileEntry &);

        static bool LessThan(const FileEntry &f1, const FileEntry &f2);

        std::string   name;
        int           size;

        bool          CanAccess() const    { return (flag & 1) == 1; }
        void          SetAccess(bool val)  { flag = (flag & 0xfe | (val?1:0)); }
        bool          IsVirtual() const    { return (flag & 2) == 2; }
        void          SetVirtual(bool val) { flag = (flag & 0xfd | (val?2:0)); }
    private:
        // Access  is bit 0
        // Virtual is bit 1
        unsigned char flag;
    };

    typedef std::vector<FileEntry> FileEntryVector;

    struct MDSERVER_PROXY_API FileList
    {
        FileList();
        FileList(const FileList &);
        ~FileList();
        void operator = (const FileList &);
        void Sort();
        void Clear();

        FileEntryVector       files;
        FileEntryVector       dirs;
        FileEntryVector       others;
        StringStringVectorMap virtualFiles;
    };
public:
    MDServerProxy();
    virtual ~MDServerProxy();

    virtual std::string GetComponentName() const;

    // RPCs to access functionality on the mdserver.
    void                       ChangeDirectory(const std::string &);
    void                       Connect(const std::string &host,
                                       const std::string &key, int port,
                                       int numRead, int numWrite);
    void                       CreateGroupList(const std::string &filename,
                                               const stringVector &groupList);
    std::string                GetDirectory();
    const FileList            *GetFileList(const std::string &filter, bool);
    const avtDatabaseMetaData *GetMetaData(const std::string &, int=0);
    const SILAttributes       *GetSIL(const std::string &, int=0);
    std::string                ExpandPath(const std::string &);
    void                       CloseDatabase();

    char                       GetSeparator() const;
    std::string                GetSeparatorString() const;
protected:
    virtual void               SetupComponentRPCs();
private:
    void                       DetermineSeparator();

    char                       separator;

    GetDirectoryRPC            getDirectoryRPC;
    ChangeDirectoryRPC         changeDirectoryRPC;
    GetFileListRPC             getFileListRPC;
    GetMetaDataRPC             getMetaDataRPC;
    GetSILRPC                  getSILRPC;
    ConnectRPC                 connectRPC;
    CreateGroupListRPC         createGroupListRPC;
    ExpandPathRPC              expandPathRPC;
    CloseDatabaseRPC           closeDatabaseRPC;

    FileList                   fileList;
    avtDatabaseMetaData        metaData;
    SILAttributes              sil;
};

#endif
