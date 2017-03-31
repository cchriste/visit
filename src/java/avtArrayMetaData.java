// ***************************************************************************
//
// Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-442911
// All rights reserved.
//
// This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
// full copyright notice is contained in the file COPYRIGHT located at the root
// of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
//
// Redistribution  and  use  in  source  and  binary  forms,  with  or  without
// modification, are permitted provided that the following conditions are met:
//
//  - Redistributions of  source code must  retain the above  copyright notice,
//    this list of conditions and the disclaimer below.
//  - Redistributions in binary form must reproduce the above copyright notice,
//    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
//    documentation and/or other materials provided with the distribution.
//  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
//    be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
// LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
// DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
// SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
// CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
// LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
// OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ***************************************************************************

package llnl.visit;

import java.util.Vector;

// ****************************************************************************
// Class: avtArrayMetaData
//
// Purpose:
//    Contains array metadata attributes
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class avtArrayMetaData extends avtVarMetaData
{
    private static int avtArrayMetaData_numAdditionalAtts = 2;

    public avtArrayMetaData()
    {
        super(avtArrayMetaData_numAdditionalAtts);

        nVars = 0;
        compNames = new Vector();
    }

    public avtArrayMetaData(int nMoreFields)
    {
        super(avtArrayMetaData_numAdditionalAtts + nMoreFields);

        nVars = 0;
        compNames = new Vector();
    }

    public avtArrayMetaData(avtArrayMetaData obj)
    {
        super(obj);

        int i;

        nVars = obj.nVars;
        compNames = new Vector(obj.compNames.size());
        for(i = 0; i < obj.compNames.size(); ++i)
            compNames.addElement(new String((String)obj.compNames.elementAt(i)));


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return avtArrayMetaData_numAdditionalAtts;
    }

    public boolean equals(avtArrayMetaData obj)
    {
        int i;

        // Compare the elements in the compNames vector.
        boolean compNames_equal = (obj.compNames.size() == compNames.size());
        for(i = 0; (i < compNames.size()) && compNames_equal; ++i)
        {
            // Make references to String from Object.
            String compNames1 = (String)compNames.elementAt(i);
            String compNames2 = (String)obj.compNames.elementAt(i);
            compNames_equal = compNames1.equals(compNames2);
        }
        // Create the return value
        return (super.equals(obj) && (nVars == obj.nVars) &&
                compNames_equal);
    }

    // Property setting methods
    public void SetNVars(int nVars_)
    {
        nVars = nVars_;
        Select((new avtArrayMetaData()).Offset() + 0);
    }

    public void SetCompNames(Vector compNames_)
    {
        compNames = compNames_;
        Select((new avtArrayMetaData()).Offset() + 1);
    }

    // Property getting methods
    public int    GetNVars() { return nVars; }
    public Vector GetCompNames() { return compNames; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        super.WriteAtts(buf);

        int offset = (new avtArrayMetaData()).Offset();
        if(WriteSelect(offset + 0, buf))
            buf.WriteInt(nVars);
        if(WriteSelect(offset + 1, buf))
            buf.WriteStringVector(compNames);
    }

    public void ReadAtts(int id, CommunicationBuffer buf)
    {
        int offset = (new avtArrayMetaData()).Offset();
        int index = id - offset;
        switch(index)
        {
        case 0:
            SetNVars(buf.ReadInt());
            break;
        case 1:
            SetCompNames(buf.ReadStringVector());
            break;
        default:
            super.ReadAtts(id, buf);
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intToString("nVars", nVars, indent) + "\n";
        str = str + stringVectorToString("compNames", compNames, indent) + "\n";
        return super.toString(indent) + str;
    }


    // Attributes
    private int    nVars;
    private Vector compNames; // vector of String objects
}

