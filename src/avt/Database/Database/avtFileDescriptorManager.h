/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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

// ************************************************************************* //
//                         avtFileDescriptorManager.h                        //
// ************************************************************************* //

#ifndef AVT_FILE_DESCRIPTOR_MANAGER_H
#define AVT_FILE_DESCRIPTOR_MANAGER_H


typedef void   (*CloseFileCallback)(void *, int);

#include <database_exports.h>

#include <vector>


// ****************************************************************************
//  Class: avtFileDescriptorManager
//
//  Purpose:
//      Manages all the open file descriptors.  This is a mechanism made use
//      of by the database and exists so that different databases don't need
//      to know about each other.
//
//  Programmer: Hank Childs
//  Creation:   March 21, 2002
//
// ****************************************************************************

class DATABASE_API avtFileDescriptorManager
{
  public:
    static avtFileDescriptorManager   *Instance(void);
    static void                        DeleteInstance(void);

    int                                RegisterFile(CloseFileCallback, void *);
    void                               UnregisterFile(int);
    void                               UsedFile(int);

    void                               SetMaximumNumberOfOpenFiles(int);

  protected:
                                       avtFileDescriptorManager();
    virtual                           ~avtFileDescriptorManager();

    int                                maximumNumberOfOpenFiles;
    int                                currentNumberOfOpenFiles;
    int                                timestamp;
    static avtFileDescriptorManager   *instance;

    std::vector<void *>                closeFileArgs;
    std::vector<CloseFileCallback>     closeFileCallbacks;
    std::vector<bool>                  fileIsOpen;
    std::vector<int>                   fileTimestamp;

    void                               CloseLeastRecentlyUsedFile(void);
};

#endif

