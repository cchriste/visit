#ifndef _GETFILELIST_RPC_H_
#define _GETFILELIST_RPC_H_
#include <mdsrpc_exports.h>

#include <VisItRPC.h>
#include <vector>
#include <string>

// ****************************************************************************
// Class: GetFileListRPC
//
// Purpose:
//   This class encapsulates a call to get a list of files from a
//   remote file system.
//
// Notes:      
//
// Programmer: Jeremy Meredith
// Creation:   Tue Aug 29 10:45:58 PDT 2000
//
// Modifications:
//   Brad Whitlock, Tue Aug 29 10:46:53 PDT 2000
//   I moved the definitions of FileList's methods to the .C file.
//
//   Brad Whitlock, Mon Mar 24 14:15:14 PST 2003
//   I added a filter string and a boolean flag to the invokation method.
//   I also added a new VIRTUAL file type. I also modified FileList so that
//   it has methods to clear and sort itself.
//
// ****************************************************************************

class MDSERVER_RPC_API GetFileListRPC : public BlockingRPC
{
public:
    enum file_types
    {
        DIR,
        REG,
        VIRTUAL,
        UNKNOWN,
    };

    struct MDSERVER_RPC_API FileList : public AttributeSubject
    {
        stringVector names;
        intVector    types;
        longVector   sizes;
        intVector    access;
        stringVector virtualNames;
        intVector    numVirtualFiles;
    public:
        FileList();
        virtual ~FileList();
        virtual void SelectAll();
        virtual const std::string TypeName() const
           { return "GetFileListRPC::FileList"; };

        void Clear();
        void Sort();
    };
public:
    GetFileListRPC();
    virtual ~GetFileListRPC();

    // Invokation method
    const FileList *operator()(const std::string &f, bool);

    // Property selection methods
    virtual void SelectAll();

    const std::string &GetFilter() const;
    bool               GetAutomaticFileGrouping() const;
private:
    FileList    fileList;
    std::string filter;
    bool        automaticFileGrouping;
};

// Method to print the file list.
MDSERVER_RPC_API ostream &operator << (ostream &os,
                                       const GetFileListRPC::FileList &fl);

#endif
