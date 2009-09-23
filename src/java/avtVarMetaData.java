// ***************************************************************************
//
// Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
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

package llnl.visit;


// ****************************************************************************
// Class: avtVarMetaData
//
// Purpose:
//    Contains metadata attributes associated with all mesh variables
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class avtVarMetaData extends avtBaseVarMetaData
{
    private static int numAdditionalAttributes = 6;

    public avtVarMetaData()
    {
        super(numAdditionalAttributes);

        centering = 0;
        hasUnits = false;
        units = new String("");
        hasDataExtents = false;
        minDataExtents = 0;
        maxDataExtents = 0;
    }

    public avtVarMetaData(int nMoreFields)
    {
        super(numAdditionalAttributes + nMoreFields);

        centering = 0;
        hasUnits = false;
        units = new String("");
        hasDataExtents = false;
        minDataExtents = 0;
        maxDataExtents = 0;
    }

    public avtVarMetaData(avtVarMetaData obj)
    {
        super(numAdditionalAttributes);

        centering = obj.centering;
        hasUnits = obj.hasUnits;
        units = new String(obj.units);
        hasDataExtents = obj.hasDataExtents;
        minDataExtents = obj.minDataExtents;
        maxDataExtents = obj.maxDataExtents;

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

    public boolean equals(avtVarMetaData obj)
    {
        // Create the return value
        return (super.equals(obj) && (centering == obj.centering) &&
                (hasUnits == obj.hasUnits) &&
                (units.equals(obj.units)) &&
                (hasDataExtents == obj.hasDataExtents) &&
                (minDataExtents == obj.minDataExtents) &&
                (maxDataExtents == obj.maxDataExtents));
    }

    // Property setting methods
    public void SetCentering(int centering_)
    {
        centering = centering_;
        Select(Offset() + 0);
    }

    public void SetHasUnits(boolean hasUnits_)
    {
        hasUnits = hasUnits_;
        Select(Offset() + 1);
    }

    public void SetUnits(String units_)
    {
        units = units_;
        Select(Offset() + 2);
    }

    public void SetHasDataExtents(boolean hasDataExtents_)
    {
        hasDataExtents = hasDataExtents_;
        Select(Offset() + 3);
    }

    public void SetMinDataExtents(double minDataExtents_)
    {
        minDataExtents = minDataExtents_;
        Select(Offset() + 4);
    }

    public void SetMaxDataExtents(double maxDataExtents_)
    {
        maxDataExtents = maxDataExtents_;
        Select(Offset() + 5);
    }

    // Property getting methods
    public int     GetCentering() { return centering; }
    public boolean GetHasUnits() { return hasUnits; }
    public String  GetUnits() { return units; }
    public boolean GetHasDataExtents() { return hasDataExtents; }
    public double  GetMinDataExtents() { return minDataExtents; }
    public double  GetMaxDataExtents() { return maxDataExtents; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        super.WriteAtts(buf);

        int offset = Offset();
        if(WriteSelect(offset + 0, buf))
            buf.WriteInt(centering);
        if(WriteSelect(offset + 1, buf))
            buf.WriteBool(hasUnits);
        if(WriteSelect(offset + 2, buf))
            buf.WriteString(units);
        if(WriteSelect(offset + 3, buf))
            buf.WriteBool(hasDataExtents);
        if(WriteSelect(offset + 4, buf))
            buf.WriteDouble(minDataExtents);
        if(WriteSelect(offset + 5, buf))
            buf.WriteDouble(maxDataExtents);
    }

    public void ReadAtts(int id, CommunicationBuffer buf)
    {
        int index = id - Offset();
        switch(index)
        {
        case 0:
            SetCentering(buf.ReadInt());
            break;
        case 1:
            SetHasUnits(buf.ReadBool());
            break;
        case 2:
            SetUnits(buf.ReadString());
            break;
        case 3:
            SetHasDataExtents(buf.ReadBool());
            break;
        case 4:
            SetMinDataExtents(buf.ReadDouble());
            break;
        case 5:
            SetMaxDataExtents(buf.ReadDouble());
            break;
        default:
            super.ReadAtts(id, buf);
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intToString("centering", centering, indent) + "\n";
        str = str + boolToString("hasUnits", hasUnits, indent) + "\n";
        str = str + stringToString("units", units, indent) + "\n";
        str = str + boolToString("hasDataExtents", hasDataExtents, indent) + "\n";
        str = str + doubleToString("minDataExtents", minDataExtents, indent) + "\n";
        str = str + doubleToString("maxDataExtents", maxDataExtents, indent) + "\n";
        return super.toString(indent) + str;
    }


    // Attributes
    private int     centering;
    private boolean hasUnits;
    private String  units;
    private boolean hasDataExtents;
    private double  minDataExtents;
    private double  maxDataExtents;
}

