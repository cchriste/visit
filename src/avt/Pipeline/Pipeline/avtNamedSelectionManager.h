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

#ifndef AVT_NAMED_SELECTION_MANAGER_H
#define AVT_NAMED_SELECTION_MANAGER_H

#include <pipeline_exports.h>

#include <vector>

#include <avtDataObject.h>
#include <avtNamedSelection.h>
#include <avtNamedSelectionExtension.h>

#include <MRUCache.h>
#include <SelectionProperties.h>
#include <visitstream.h>

// ****************************************************************************
//  Class: avtNamedSelectionManager
//
//  Purpose:
//      Manager for named selections.
//
//  Programmer: Hank Childs
//  Creation:   January 30, 2009
//
//  Modifications:
//
//    Hank Childs, Sun Apr 19 22:42:09 PDT 2009
//    Add argument to DeleteNamedSelection.
//
//    Hank Childs, Mon Jul 13 15:53:54 PDT 2009
//    Add arguments to [Get|Save]NamedSelection for automatic saving and
//    loading of selections to help with fault tolerance, save/restore 
//    sessions, etc.
//
//    Brad Whitlock, Mon Dec 13 16:15:41 PST 2010
//    I added an avtNamedSelectionExtension argument to CreateNamedSelection,
//    which lets us pass in an object that can perform additional setup for a
//    selection based on its properties. This class manages named selections,
//    which are the actual resulting data that is used to restrict selections.
//    I also added selection properties which describe how the selection was 
//    defined. These properties can be queried in the event that we need to
//    create a selection using some other method.
//
//    Brad Whitlock, Tue Sep  6 15:15:53 PDT 2011
//    I added a cache.
//
// ****************************************************************************

class PIPELINE_API avtNamedSelectionManager
{
  public:
                  avtNamedSelectionManager();
    virtual      ~avtNamedSelectionManager();
    
    static avtNamedSelectionManager *
                  GetInstance(void);

    avtNamedSelection *
                  GetNamedSelection(const std::string &);

    void          CreateNamedSelection(avtDataObject_p,
                                       const SelectionProperties &,
                                       avtNamedSelectionExtension *);

    void          DeleteNamedSelection(const std::string &, 
                                       bool expectThisSelToBeThere);
    bool          LoadNamedSelection(const std::string &, bool = false);
    void          SaveNamedSelection(const std::string &, bool = false);

    void          ClearCache(const std::string &selName = std::string());

    const SelectionProperties *GetSelectionProperties(const std::string &selName) const;

    static int MaximumSelectionSize();

  protected:
    static avtNamedSelectionManager    *instance;
    std::vector<avtNamedSelection *>    selList;

    void          AddSelectionProperties(const SelectionProperties &);
    std::vector<SelectionProperties>    properties;

    avtNamedSelectionCache              cache;

  private:
    // These methods are defined to prevent accidental use of bitwise copy
    // implementations.  If you want to re-define them to do something
    // meaningful, that's fine.
                         avtNamedSelectionManager(const avtNamedSelectionManager &) {;};
    avtNamedSelectionManager          &operator=(const avtNamedSelectionManager &) { return *this; };

    std::string          CreateQualifiedSelectionName(const std::string &, bool);
    avtNamedSelection   *IterateOverNamedSelections(const std::string &);
};


#endif


