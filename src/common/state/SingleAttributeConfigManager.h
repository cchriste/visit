/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
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

#ifndef SINGLE_ATTRIBUTE_MANAGER_H
#define SINGLE_ATTRIBUTE_MANAGER_H
#include <ConfigManager.h>
#include <AttributeGroup.h>
#include <string>

class DataNode;

// ****************************************************************************
// Class: SingleAttributeConfigManager
//
// Purpose:
//   Reads and writes single attributes to/from files.
//
// Notes:      
//
// Programmer: Jeremy Meredith
// Creation:   January  2, 2009
//
// Modifications:
//   Jeremy Meredith, Thu Apr 29 12:14:13 EDT 2010
//   Added ability to do a selective save.
//
// ****************************************************************************

class STATE_API SingleAttributeConfigManager : public ConfigManager
{
public:
    SingleAttributeConfigManager(AttributeGroup*);
    virtual ~SingleAttributeConfigManager();

    bool              Export(const std::string &filename,
                             bool completeSave = true);
    bool              Import(const std::string &filename);

    bool              Export(std::ostream& out,
                             bool completeSave = true);
    bool              Import(std::istream& in);

    virtual bool      WriteConfigFile(const char *filename);
    virtual DataNode *ReadConfigFile(const char *filename);

    virtual bool      WriteConfigFile(std::ostream& out);
    virtual DataNode *ReadConfigFile(std::istream& in);
protected:
private:
    AttributeGroup *attribute;
    bool completeSave;
};

#endif
