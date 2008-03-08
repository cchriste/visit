// ***************************************************************************
//
// Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400142
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
// Class: CylinderAttributes
//
// Purpose:
//    Contain the attributes for a cylinder
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class CylinderAttributes extends AttributeSubject implements Plugin
{
    public CylinderAttributes()
    {
        super(3);

        point1 = new double[3];
        point1[0] = 0;
        point1[1] = 0;
        point1[2] = 0;
        point2 = new double[3];
        point2[0] = 1;
        point2[1] = 0;
        point2[2] = 0;
        radius = 1;
    }

    public CylinderAttributes(CylinderAttributes obj)
    {
        super(3);

        int i;

        point1 = new double[3];
        point1[0] = obj.point1[0];
        point1[1] = obj.point1[1];
        point1[2] = obj.point1[2];

        point2 = new double[3];
        point2[0] = obj.point2[0];
        point2[1] = obj.point2[1];
        point2[2] = obj.point2[2];

        radius = obj.radius;

        SelectAll();
    }

    public boolean equals(CylinderAttributes obj)
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
        return (point1_equal &&
                point2_equal &&
                (radius == obj.radius));
    }

    public String GetName() { return "Cylinder"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetPoint1(double[] point1_)
    {
        point1[0] = point1_[0];
        point1[1] = point1_[1];
        point1[2] = point1_[2];
        Select(0);
    }

    public void SetPoint1(double e0, double e1, double e2)
    {
        point1[0] = e0;
        point1[1] = e1;
        point1[2] = e2;
        Select(0);
    }

    public void SetPoint2(double[] point2_)
    {
        point2[0] = point2_[0];
        point2[1] = point2_[1];
        point2[2] = point2_[2];
        Select(1);
    }

    public void SetPoint2(double e0, double e1, double e2)
    {
        point2[0] = e0;
        point2[1] = e1;
        point2[2] = e2;
        Select(1);
    }

    public void SetRadius(double radius_)
    {
        radius = radius_;
        Select(2);
    }

    // Property getting methods
    public double[] GetPoint1() { return point1; }
    public double[] GetPoint2() { return point2; }
    public double   GetRadius() { return radius; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(point1);
        if(WriteSelect(1, buf))
            buf.WriteDoubleArray(point2);
        if(WriteSelect(2, buf))
            buf.WriteDouble(radius);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetPoint1(buf.ReadDoubleArray());
                break;
            case 1:
                SetPoint2(buf.ReadDoubleArray());
                break;
            case 2:
                SetRadius(buf.ReadDouble());
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleArrayToString("point1", point1, indent) + "\n";
        str = str + doubleArrayToString("point2", point2, indent) + "\n";
        str = str + doubleToString("radius", radius, indent) + "\n";
        return str;
    }


    // Attributes
    private double[] point1;
    private double[] point2;
    private double   radius;
}

