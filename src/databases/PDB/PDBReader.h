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

#ifndef PDB_READER_H
#define PDB_READER_H
#include <string>
#include <map>
#include <PDBFileObject.h>

// ****************************************************************************
// Class: PDBReader
//
// Purpose:
//   Abstract base class for PDB file format readers.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Tue Sep 16 09:25:28 PDT 2003
//
// Modifications:
//   Brad Whitlock, Thu Sep 2 00:06:54 PDT 2004
//   Added the VariableData::FreeData method.
//
//   Brad Whitlock, Thu Nov  6 14:39:36 PST 2008
//   Added method to get the PDB file.
//
//   Mark C. Miller, Tue Apr 28 11:05:54 PDT 2009
//   Changed name of PDB() to PDBfobj() to avoid symbol collision with PDB
//   proper.
// ****************************************************************************

class PDBReader
{
public:
    PDBReader(const char *filename);
    PDBReader(PDBFileObject *p);
    virtual ~PDBReader();
    void Close();
    void SetOwnsPDBFile(bool v);

    bool Identify();
    PDBFileObject *PDBfobj();
protected:
    class VariableData
    {
    public:
        VariableData(const std::string &name);
        ~VariableData();
        bool ReadValues(PDBFileObject *pdb);
        void FreeData();

        std::string  varName;
        void        *data;
        TypeEnum     dataType;
        int         *dims;
        int          nDims;
        int          nTotalElements;
    };

    typedef std::map<std::string, VariableData *> VariableDataMap;

    virtual bool IdentifyFormat() = 0;

    PDBFileObject *pdb;
    bool           ownsPDBFile;
};

#endif
