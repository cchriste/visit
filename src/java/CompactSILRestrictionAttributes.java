// ***************************************************************************
//
// Copyright (c) 2000 - 2012, Lawrence Livermore National Security, LLC
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

import java.lang.Byte;
import java.util.Vector;

// ****************************************************************************
// Class: CompactSILRestrictionAttributes
//
// Purpose:
//    The class contains attributes for a compacted SIL restrictions.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class CompactSILRestrictionAttributes extends AttributeSubject
{
    private static int CompactSILRestrictionAttributes_numAdditionalAtts = 3;

    public CompactSILRestrictionAttributes()
    {
        super(CompactSILRestrictionAttributes_numAdditionalAtts);

        useSet = new Vector();
        topSet = new String("");
        topSetIsAllOn = false;
    }

    public CompactSILRestrictionAttributes(int nMoreFields)
    {
        super(CompactSILRestrictionAttributes_numAdditionalAtts + nMoreFields);

        useSet = new Vector();
        topSet = new String("");
        topSetIsAllOn = false;
    }

    public CompactSILRestrictionAttributes(CompactSILRestrictionAttributes obj)
    {
        super(CompactSILRestrictionAttributes_numAdditionalAtts);

        int i;

        useSet = new Vector(obj.useSet.size());
        for(i = 0; i < obj.useSet.size(); ++i)
        {
            Byte bv = (Byte)obj.useSet.elementAt(i);
            useSet.addElement(new Byte(bv.byteValue()));
        }

        topSet = new String(obj.topSet);
        topSetIsAllOn = obj.topSetIsAllOn;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return CompactSILRestrictionAttributes_numAdditionalAtts;
    }

    public boolean equals(CompactSILRestrictionAttributes obj)
    {
        int i;

        // Compare the elements in the useSet vector.
        boolean useSet_equal = (obj.useSet.size() == useSet.size());
        for(i = 0; (i < useSet.size()) && useSet_equal; ++i)
        {
            // Make references to Byte from Object.
            Byte useSet1 = (Byte)useSet.elementAt(i);
            Byte useSet2 = (Byte)obj.useSet.elementAt(i);
            useSet_equal = useSet1.equals(useSet2);
        }
        // Create the return value
        return (useSet_equal &&
                (topSet.equals(obj.topSet)) &&
                (topSetIsAllOn == obj.topSetIsAllOn));
    }

    // Property setting methods
    public void SetUseSet(Vector useSet_)
    {
        useSet = useSet_;
        Select(0);
    }

    public void SetTopSet(String topSet_)
    {
        topSet = topSet_;
        Select(1);
    }

    public void SetTopSetIsAllOn(boolean topSetIsAllOn_)
    {
        topSetIsAllOn = topSetIsAllOn_;
        Select(2);
    }

    // Property getting methods
    public Vector  GetUseSet() { return useSet; }
    public String  GetTopSet() { return topSet; }
    public boolean GetTopSetIsAllOn() { return topSetIsAllOn; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteByteVector(useSet);
        if(WriteSelect(1, buf))
            buf.WriteString(topSet);
        if(WriteSelect(2, buf))
            buf.WriteBool(topSetIsAllOn);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetUseSet(buf.ReadByteVector());
            break;
        case 1:
            SetTopSet(buf.ReadString());
            break;
        case 2:
            SetTopSetIsAllOn(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + ucharVectorToString("useSet", useSet, indent) + "\n";
        str = str + stringToString("topSet", topSet, indent) + "\n";
        str = str + boolToString("topSetIsAllOn", topSetIsAllOn, indent) + "\n";
        return str;
    }


    // Attributes
    private Vector  useSet; // vector of Byte objects
    private String  topSet;
    private boolean topSetIsAllOn;
}

