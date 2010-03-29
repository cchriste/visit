// ***************************************************************************
//
// Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400124
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

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: CoordSwapAttributes
//
// Purpose:
//    This class contains attributes for the coord swap operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class CoordSwapAttributes extends AttributeSubject implements Plugin
{
    private static int numAdditionalAttributes = 3;

    // Enum values
    public final static int COORD_COORD1 = 0;
    public final static int COORD_COORD2 = 1;
    public final static int COORD_COORD3 = 2;


    public CoordSwapAttributes()
    {
        super(numAdditionalAttributes);

        newCoord1 = COORD_COORD1;
        newCoord2 = COORD_COORD2;
        newCoord3 = COORD_COORD3;
    }

    public CoordSwapAttributes(int nMoreFields)
    {
        super(numAdditionalAttributes + nMoreFields);

        newCoord1 = COORD_COORD1;
        newCoord2 = COORD_COORD2;
        newCoord3 = COORD_COORD3;
    }

    public CoordSwapAttributes(CoordSwapAttributes obj)
    {
        super(numAdditionalAttributes);

        newCoord1 = obj.newCoord1;
        newCoord2 = obj.newCoord2;
        newCoord3 = obj.newCoord3;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return numAdditionalAttributes;
    }

    public boolean equals(CoordSwapAttributes obj)
    {
        // Create the return value
        return ((newCoord1 == obj.newCoord1) &&
                (newCoord2 == obj.newCoord2) &&
                (newCoord3 == obj.newCoord3));
    }

    public String GetName() { return "CoordSwap"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetNewCoord1(int newCoord1_)
    {
        newCoord1 = newCoord1_;
        Select(0);
    }

    public void SetNewCoord2(int newCoord2_)
    {
        newCoord2 = newCoord2_;
        Select(1);
    }

    public void SetNewCoord3(int newCoord3_)
    {
        newCoord3 = newCoord3_;
        Select(2);
    }

    // Property getting methods
    public int GetNewCoord1() { return newCoord1; }
    public int GetNewCoord2() { return newCoord2; }
    public int GetNewCoord3() { return newCoord3; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(newCoord1);
        if(WriteSelect(1, buf))
            buf.WriteInt(newCoord2);
        if(WriteSelect(2, buf))
            buf.WriteInt(newCoord3);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetNewCoord1(buf.ReadInt());
            break;
        case 1:
            SetNewCoord2(buf.ReadInt());
            break;
        case 2:
            SetNewCoord3(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "newCoord1 = ";
        if(newCoord1 == COORD_COORD1)
            str = str + "COORD_COORD1";
        if(newCoord1 == COORD_COORD2)
            str = str + "COORD_COORD2";
        if(newCoord1 == COORD_COORD3)
            str = str + "COORD_COORD3";
        str = str + "\n";
        str = str + indent + "newCoord2 = ";
        if(newCoord2 == COORD_COORD1)
            str = str + "COORD_COORD1";
        if(newCoord2 == COORD_COORD2)
            str = str + "COORD_COORD2";
        if(newCoord2 == COORD_COORD3)
            str = str + "COORD_COORD3";
        str = str + "\n";
        str = str + indent + "newCoord3 = ";
        if(newCoord3 == COORD_COORD1)
            str = str + "COORD_COORD1";
        if(newCoord3 == COORD_COORD2)
            str = str + "COORD_COORD2";
        if(newCoord3 == COORD_COORD3)
            str = str + "COORD_COORD3";
        str = str + "\n";
        return str;
    }


    // Attributes
    private int newCoord1;
    private int newCoord2;
    private int newCoord3;
}

