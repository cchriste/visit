// ***************************************************************************
//
// Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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
// Class: LineSurfaceAttributes
//
// Purpose:
//    Attributes for the LineSurface operator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class LineSurfaceAttributes extends AttributeSubject implements Plugin
{
    private static int LineSurfaceAttributes_numAdditionalAtts = 5;

    public LineSurfaceAttributes()
    {
        super(LineSurfaceAttributes_numAdditionalAtts);

        startTime = 0;
        endTime = 1;
        stride = 1;
        point1 = new double[3];
        point1[0] = 0;
        point1[1] = 0;
        point1[2] = 0;
        point2 = new double[3];
        point2[0] = 0;
        point2[1] = 0;
        point2[2] = 0;
    }

    public LineSurfaceAttributes(int nMoreFields)
    {
        super(LineSurfaceAttributes_numAdditionalAtts + nMoreFields);

        startTime = 0;
        endTime = 1;
        stride = 1;
        point1 = new double[3];
        point1[0] = 0;
        point1[1] = 0;
        point1[2] = 0;
        point2 = new double[3];
        point2[0] = 0;
        point2[1] = 0;
        point2[2] = 0;
    }

    public LineSurfaceAttributes(LineSurfaceAttributes obj)
    {
        super(LineSurfaceAttributes_numAdditionalAtts);

        int i;

        startTime = obj.startTime;
        endTime = obj.endTime;
        stride = obj.stride;
        point1 = new double[3];
        point1[0] = obj.point1[0];
        point1[1] = obj.point1[1];
        point1[2] = obj.point1[2];

        point2 = new double[3];
        point2[0] = obj.point2[0];
        point2[1] = obj.point2[1];
        point2[2] = obj.point2[2];


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return LineSurfaceAttributes_numAdditionalAtts;
    }

    public boolean equals(LineSurfaceAttributes obj)
    {
        int i;

        // Compare the point1 arrays.
        boolean point1_equal = true;
        for(i = 0; i < 3 && point1_equal; ++i)
            point1_equal = (point1[i] == obj.point1[i]);

        // Compare the point2 arrays.
        boolean point2_equal = true;
        for(i = 0; i < 3 && point2_equal; ++i)
            point2_equal = (point2[i] == obj.point2[i]);

        // Create the return value
        return ((startTime == obj.startTime) &&
                (endTime == obj.endTime) &&
                (stride == obj.stride) &&
                point1_equal &&
                point2_equal);
    }

    public String GetName() { return "LineSurface"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetStartTime(int startTime_)
    {
        startTime = startTime_;
        Select(0);
    }

    public void SetEndTime(int endTime_)
    {
        endTime = endTime_;
        Select(1);
    }

    public void SetStride(int stride_)
    {
        stride = stride_;
        Select(2);
    }

    public void SetPoint1(double[] point1_)
    {
        point1[0] = point1_[0];
        point1[1] = point1_[1];
        point1[2] = point1_[2];
        Select(3);
    }

    public void SetPoint1(double e0, double e1, double e2)
    {
        point1[0] = e0;
        point1[1] = e1;
        point1[2] = e2;
        Select(3);
    }

    public void SetPoint2(double[] point2_)
    {
        point2[0] = point2_[0];
        point2[1] = point2_[1];
        point2[2] = point2_[2];
        Select(4);
    }

    public void SetPoint2(double e0, double e1, double e2)
    {
        point2[0] = e0;
        point2[1] = e1;
        point2[2] = e2;
        Select(4);
    }

    // Property getting methods
    public int      GetStartTime() { return startTime; }
    public int      GetEndTime() { return endTime; }
    public int      GetStride() { return stride; }
    public double[] GetPoint1() { return point1; }
    public double[] GetPoint2() { return point2; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(startTime);
        if(WriteSelect(1, buf))
            buf.WriteInt(endTime);
        if(WriteSelect(2, buf))
            buf.WriteInt(stride);
        if(WriteSelect(3, buf))
            buf.WriteDoubleArray(point1);
        if(WriteSelect(4, buf))
            buf.WriteDoubleArray(point2);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetStartTime(buf.ReadInt());
            break;
        case 1:
            SetEndTime(buf.ReadInt());
            break;
        case 2:
            SetStride(buf.ReadInt());
            break;
        case 3:
            SetPoint1(buf.ReadDoubleArray());
            break;
        case 4:
            SetPoint2(buf.ReadDoubleArray());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intToString("startTime", startTime, indent) + "\n";
        str = str + intToString("endTime", endTime, indent) + "\n";
        str = str + intToString("stride", stride, indent) + "\n";
        str = str + doubleArrayToString("point1", point1, indent) + "\n";
        str = str + doubleArrayToString("point2", point2, indent) + "\n";
        return str;
    }


    // Attributes
    private int      startTime;
    private int      endTime;
    private int      stride;
    private double[] point1;
    private double[] point2;
}

