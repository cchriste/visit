// ***************************************************************************
//
// Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: ExtrudeAttributes
//
// Purpose:
//    This class contains attributes for the extrude operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ExtrudeAttributes extends AttributeSubject implements Plugin
{
    private static int ExtrudeAttributes_numAdditionalAtts = 4;

    public ExtrudeAttributes()
    {
        super(ExtrudeAttributes_numAdditionalAtts);

        axis = new double[3];
        axis[0] = 0;
        axis[1] = 0;
        axis[2] = 1;
        length = 1;
        steps = 30;
        preserveOriginalCellNumbers = true;
    }

    public ExtrudeAttributes(int nMoreFields)
    {
        super(ExtrudeAttributes_numAdditionalAtts + nMoreFields);

        axis = new double[3];
        axis[0] = 0;
        axis[1] = 0;
        axis[2] = 1;
        length = 1;
        steps = 30;
        preserveOriginalCellNumbers = true;
    }

    public ExtrudeAttributes(ExtrudeAttributes obj)
    {
        super(ExtrudeAttributes_numAdditionalAtts);

        int i;

        axis = new double[3];
        axis[0] = obj.axis[0];
        axis[1] = obj.axis[1];
        axis[2] = obj.axis[2];

        length = obj.length;
        steps = obj.steps;
        preserveOriginalCellNumbers = obj.preserveOriginalCellNumbers;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ExtrudeAttributes_numAdditionalAtts;
    }

    public boolean equals(ExtrudeAttributes obj)
    {
        int i;

        // Compare the axis arrays.
        boolean axis_equal = true;
        for(i = 0; i < 3 && axis_equal; ++i)
            axis_equal = (axis[i] == obj.axis[i]);

        // Create the return value
        return (axis_equal &&
                (length == obj.length) &&
                (steps == obj.steps) &&
                (preserveOriginalCellNumbers == obj.preserveOriginalCellNumbers));
    }

    public String GetName() { return "Extrude"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetAxis(double[] axis_)
    {
        axis[0] = axis_[0];
        axis[1] = axis_[1];
        axis[2] = axis_[2];
        Select(0);
    }

    public void SetAxis(double e0, double e1, double e2)
    {
        axis[0] = e0;
        axis[1] = e1;
        axis[2] = e2;
        Select(0);
    }

    public void SetLength(double length_)
    {
        length = length_;
        Select(1);
    }

    public void SetSteps(int steps_)
    {
        steps = steps_;
        Select(2);
    }

    public void SetPreserveOriginalCellNumbers(boolean preserveOriginalCellNumbers_)
    {
        preserveOriginalCellNumbers = preserveOriginalCellNumbers_;
        Select(3);
    }

    // Property getting methods
    public double[] GetAxis() { return axis; }
    public double   GetLength() { return length; }
    public int      GetSteps() { return steps; }
    public boolean  GetPreserveOriginalCellNumbers() { return preserveOriginalCellNumbers; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(axis);
        if(WriteSelect(1, buf))
            buf.WriteDouble(length);
        if(WriteSelect(2, buf))
            buf.WriteInt(steps);
        if(WriteSelect(3, buf))
            buf.WriteBool(preserveOriginalCellNumbers);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetAxis(buf.ReadDoubleArray());
            break;
        case 1:
            SetLength(buf.ReadDouble());
            break;
        case 2:
            SetSteps(buf.ReadInt());
            break;
        case 3:
            SetPreserveOriginalCellNumbers(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleArrayToString("axis", axis, indent) + "\n";
        str = str + doubleToString("length", length, indent) + "\n";
        str = str + intToString("steps", steps, indent) + "\n";
        str = str + boolToString("preserveOriginalCellNumbers", preserveOriginalCellNumbers, indent) + "\n";
        return str;
    }


    // Attributes
    private double[] axis;
    private double   length;
    private int      steps;
    private boolean  preserveOriginalCellNumbers;
}

