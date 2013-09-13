// ***************************************************************************
//
// Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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


// ****************************************************************************
// Class: AxisTickMarks
//
// Purpose:
//    Contains the tick mark properties for one axis.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class AxisTickMarks extends AttributeSubject
{
    private static int AxisTickMarks_numAdditionalAtts = 5;

    public AxisTickMarks()
    {
        super(AxisTickMarks_numAdditionalAtts);

        visible = true;
        majorMinimum = 0;
        majorMaximum = 1;
        minorSpacing = 0.02;
        majorSpacing = 0.2;
    }

    public AxisTickMarks(int nMoreFields)
    {
        super(AxisTickMarks_numAdditionalAtts + nMoreFields);

        visible = true;
        majorMinimum = 0;
        majorMaximum = 1;
        minorSpacing = 0.02;
        majorSpacing = 0.2;
    }

    public AxisTickMarks(AxisTickMarks obj)
    {
        super(AxisTickMarks_numAdditionalAtts);

        visible = obj.visible;
        majorMinimum = obj.majorMinimum;
        majorMaximum = obj.majorMaximum;
        minorSpacing = obj.minorSpacing;
        majorSpacing = obj.majorSpacing;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return AxisTickMarks_numAdditionalAtts;
    }

    public boolean equals(AxisTickMarks obj)
    {
        // Create the return value
        return ((visible == obj.visible) &&
                (majorMinimum == obj.majorMinimum) &&
                (majorMaximum == obj.majorMaximum) &&
                (minorSpacing == obj.minorSpacing) &&
                (majorSpacing == obj.majorSpacing));
    }

    // Property setting methods
    public void SetVisible(boolean visible_)
    {
        visible = visible_;
        Select(0);
    }

    public void SetMajorMinimum(double majorMinimum_)
    {
        majorMinimum = majorMinimum_;
        Select(1);
    }

    public void SetMajorMaximum(double majorMaximum_)
    {
        majorMaximum = majorMaximum_;
        Select(2);
    }

    public void SetMinorSpacing(double minorSpacing_)
    {
        minorSpacing = minorSpacing_;
        Select(3);
    }

    public void SetMajorSpacing(double majorSpacing_)
    {
        majorSpacing = majorSpacing_;
        Select(4);
    }

    // Property getting methods
    public boolean GetVisible() { return visible; }
    public double  GetMajorMinimum() { return majorMinimum; }
    public double  GetMajorMaximum() { return majorMaximum; }
    public double  GetMinorSpacing() { return minorSpacing; }
    public double  GetMajorSpacing() { return majorSpacing; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(visible);
        if(WriteSelect(1, buf))
            buf.WriteDouble(majorMinimum);
        if(WriteSelect(2, buf))
            buf.WriteDouble(majorMaximum);
        if(WriteSelect(3, buf))
            buf.WriteDouble(minorSpacing);
        if(WriteSelect(4, buf))
            buf.WriteDouble(majorSpacing);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetVisible(buf.ReadBool());
            break;
        case 1:
            SetMajorMinimum(buf.ReadDouble());
            break;
        case 2:
            SetMajorMaximum(buf.ReadDouble());
            break;
        case 3:
            SetMinorSpacing(buf.ReadDouble());
            break;
        case 4:
            SetMajorSpacing(buf.ReadDouble());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + boolToString("visible", visible, indent) + "\n";
        str = str + doubleToString("majorMinimum", majorMinimum, indent) + "\n";
        str = str + doubleToString("majorMaximum", majorMaximum, indent) + "\n";
        str = str + doubleToString("minorSpacing", minorSpacing, indent) + "\n";
        str = str + doubleToString("majorSpacing", majorSpacing, indent) + "\n";
        return str;
    }


    // Attributes
    private boolean visible;
    private double  majorMinimum;
    private double  majorMaximum;
    private double  minorSpacing;
    private double  majorSpacing;
}

