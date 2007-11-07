// ***************************************************************************
//
// Copyright (c) 2000 - 2007, The Regents of the University of California
// Produced at the Lawrence Livermore National Laboratory
// All rights reserved.
//
// This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
//    documentation and/or materials provided with the distribution.
//  - Neither the name of the UC/LLNL nor  the names of its contributors may be
//    used to  endorse or  promote products derived from  this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
// CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
// ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
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
// Class: ClipAttributes
//
// Purpose:
//    This class contains attributes for the clip operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Nov 7 14:50:46 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class ClipAttributes extends AttributeSubject implements Plugin
{
    // Enum values
    public final static int CLIPSTYLE_PLANE = 0;
    public final static int CLIPSTYLE_SPHERE = 1;

    public final static int WHICHCLIPPLANE_NONE = 0;
    public final static int WHICHCLIPPLANE_PLANE1 = 1;
    public final static int WHICHCLIPPLANE_PLANE2 = 2;
    public final static int WHICHCLIPPLANE_PLANE3 = 3;


    public ClipAttributes()
    {
        super(15);

        funcType = CLIPSTYLE_PLANE;
        plane1Status = true;
        plane2Status = false;
        plane3Status = false;
        plane1Origin = new double[3];
        plane1Origin[0] = 0;
        plane1Origin[1] = 0;
        plane1Origin[2] = 0;
        plane2Origin = new double[3];
        plane2Origin[0] = 0;
        plane2Origin[1] = 0;
        plane2Origin[2] = 0;
        plane3Origin = new double[3];
        plane3Origin[0] = 0;
        plane3Origin[1] = 0;
        plane3Origin[2] = 0;
        plane1Normal = new double[3];
        plane1Normal[0] = 1;
        plane1Normal[1] = 0;
        plane1Normal[2] = 0;
        plane2Normal = new double[3];
        plane2Normal[0] = 0;
        plane2Normal[1] = 1;
        plane2Normal[2] = 0;
        plane3Normal = new double[3];
        plane3Normal[0] = 0;
        plane3Normal[1] = 0;
        plane3Normal[2] = 1;
        planeInverse = false;
        planeToolControlledClipPlane = WHICHCLIPPLANE_PLANE1;
        center = new double[3];
        center[0] = 0;
        center[1] = 0;
        center[2] = 0;
        radius = 1;
        sphereInverse = false;
    }

    public ClipAttributes(ClipAttributes obj)
    {
        super(15);

        int i;

        funcType = obj.funcType;
        plane1Status = obj.plane1Status;
        plane2Status = obj.plane2Status;
        plane3Status = obj.plane3Status;
        plane1Origin = new double[3];
        plane1Origin[0] = obj.plane1Origin[0];
        plane1Origin[1] = obj.plane1Origin[1];
        plane1Origin[2] = obj.plane1Origin[2];

        plane2Origin = new double[3];
        plane2Origin[0] = obj.plane2Origin[0];
        plane2Origin[1] = obj.plane2Origin[1];
        plane2Origin[2] = obj.plane2Origin[2];

        plane3Origin = new double[3];
        plane3Origin[0] = obj.plane3Origin[0];
        plane3Origin[1] = obj.plane3Origin[1];
        plane3Origin[2] = obj.plane3Origin[2];

        plane1Normal = new double[3];
        plane1Normal[0] = obj.plane1Normal[0];
        plane1Normal[1] = obj.plane1Normal[1];
        plane1Normal[2] = obj.plane1Normal[2];

        plane2Normal = new double[3];
        plane2Normal[0] = obj.plane2Normal[0];
        plane2Normal[1] = obj.plane2Normal[1];
        plane2Normal[2] = obj.plane2Normal[2];

        plane3Normal = new double[3];
        plane3Normal[0] = obj.plane3Normal[0];
        plane3Normal[1] = obj.plane3Normal[1];
        plane3Normal[2] = obj.plane3Normal[2];

        planeInverse = obj.planeInverse;
        planeToolControlledClipPlane = obj.planeToolControlledClipPlane;
        center = new double[3];
        center[0] = obj.center[0];
        center[1] = obj.center[1];
        center[2] = obj.center[2];

        radius = obj.radius;
        sphereInverse = obj.sphereInverse;

        SelectAll();
    }

    public boolean equals(ClipAttributes obj)
    {
        int i;

        // Compare the plane1Origin arrays.
        boolean plane1Origin_equal = true;
        for(i = 0; i < 3 && plane1Origin_equal; ++i)
            plane1Origin_equal = (plane1Origin[i] == obj.plane1Origin[i]);

        // Compare the plane2Origin arrays.
        boolean plane2Origin_equal = true;
        for(i = 0; i < 3 && plane2Origin_equal; ++i)
            plane2Origin_equal = (plane2Origin[i] == obj.plane2Origin[i]);

        // Compare the plane3Origin arrays.
        boolean plane3Origin_equal = true;
        for(i = 0; i < 3 && plane3Origin_equal; ++i)
            plane3Origin_equal = (plane3Origin[i] == obj.plane3Origin[i]);

        // Compare the plane1Normal arrays.
        boolean plane1Normal_equal = true;
        for(i = 0; i < 3 && plane1Normal_equal; ++i)
            plane1Normal_equal = (plane1Normal[i] == obj.plane1Normal[i]);

        // Compare the plane2Normal arrays.
        boolean plane2Normal_equal = true;
        for(i = 0; i < 3 && plane2Normal_equal; ++i)
            plane2Normal_equal = (plane2Normal[i] == obj.plane2Normal[i]);

        // Compare the plane3Normal arrays.
        boolean plane3Normal_equal = true;
        for(i = 0; i < 3 && plane3Normal_equal; ++i)
            plane3Normal_equal = (plane3Normal[i] == obj.plane3Normal[i]);

        // Compare the center arrays.
        boolean center_equal = true;
        for(i = 0; i < 3 && center_equal; ++i)
            center_equal = (center[i] == obj.center[i]);

        // Create the return value
        return ((funcType == obj.funcType) &&
                (plane1Status == obj.plane1Status) &&
                (plane2Status == obj.plane2Status) &&
                (plane3Status == obj.plane3Status) &&
                plane1Origin_equal &&
                plane2Origin_equal &&
                plane3Origin_equal &&
                plane1Normal_equal &&
                plane2Normal_equal &&
                plane3Normal_equal &&
                (planeInverse == obj.planeInverse) &&
                (planeToolControlledClipPlane == obj.planeToolControlledClipPlane) &&
                center_equal &&
                (radius == obj.radius) &&
                (sphereInverse == obj.sphereInverse));
    }

    public String GetName() { return "Clip"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetFuncType(int funcType_)
    {
        funcType = funcType_;
        Select(0);
    }

    public void SetPlane1Status(boolean plane1Status_)
    {
        plane1Status = plane1Status_;
        Select(1);
    }

    public void SetPlane2Status(boolean plane2Status_)
    {
        plane2Status = plane2Status_;
        Select(2);
    }

    public void SetPlane3Status(boolean plane3Status_)
    {
        plane3Status = plane3Status_;
        Select(3);
    }

    public void SetPlane1Origin(double[] plane1Origin_)
    {
        plane1Origin[0] = plane1Origin_[0];
        plane1Origin[1] = plane1Origin_[1];
        plane1Origin[2] = plane1Origin_[2];
        Select(4);
    }

    public void SetPlane1Origin(double e0, double e1, double e2)
    {
        plane1Origin[0] = e0;
        plane1Origin[1] = e1;
        plane1Origin[2] = e2;
        Select(4);
    }

    public void SetPlane2Origin(double[] plane2Origin_)
    {
        plane2Origin[0] = plane2Origin_[0];
        plane2Origin[1] = plane2Origin_[1];
        plane2Origin[2] = plane2Origin_[2];
        Select(5);
    }

    public void SetPlane2Origin(double e0, double e1, double e2)
    {
        plane2Origin[0] = e0;
        plane2Origin[1] = e1;
        plane2Origin[2] = e2;
        Select(5);
    }

    public void SetPlane3Origin(double[] plane3Origin_)
    {
        plane3Origin[0] = plane3Origin_[0];
        plane3Origin[1] = plane3Origin_[1];
        plane3Origin[2] = plane3Origin_[2];
        Select(6);
    }

    public void SetPlane3Origin(double e0, double e1, double e2)
    {
        plane3Origin[0] = e0;
        plane3Origin[1] = e1;
        plane3Origin[2] = e2;
        Select(6);
    }

    public void SetPlane1Normal(double[] plane1Normal_)
    {
        plane1Normal[0] = plane1Normal_[0];
        plane1Normal[1] = plane1Normal_[1];
        plane1Normal[2] = plane1Normal_[2];
        Select(7);
    }

    public void SetPlane1Normal(double e0, double e1, double e2)
    {
        plane1Normal[0] = e0;
        plane1Normal[1] = e1;
        plane1Normal[2] = e2;
        Select(7);
    }

    public void SetPlane2Normal(double[] plane2Normal_)
    {
        plane2Normal[0] = plane2Normal_[0];
        plane2Normal[1] = plane2Normal_[1];
        plane2Normal[2] = plane2Normal_[2];
        Select(8);
    }

    public void SetPlane2Normal(double e0, double e1, double e2)
    {
        plane2Normal[0] = e0;
        plane2Normal[1] = e1;
        plane2Normal[2] = e2;
        Select(8);
    }

    public void SetPlane3Normal(double[] plane3Normal_)
    {
        plane3Normal[0] = plane3Normal_[0];
        plane3Normal[1] = plane3Normal_[1];
        plane3Normal[2] = plane3Normal_[2];
        Select(9);
    }

    public void SetPlane3Normal(double e0, double e1, double e2)
    {
        plane3Normal[0] = e0;
        plane3Normal[1] = e1;
        plane3Normal[2] = e2;
        Select(9);
    }

    public void SetPlaneInverse(boolean planeInverse_)
    {
        planeInverse = planeInverse_;
        Select(10);
    }

    public void SetPlaneToolControlledClipPlane(int planeToolControlledClipPlane_)
    {
        planeToolControlledClipPlane = planeToolControlledClipPlane_;
        Select(11);
    }

    public void SetCenter(double[] center_)
    {
        center[0] = center_[0];
        center[1] = center_[1];
        center[2] = center_[2];
        Select(12);
    }

    public void SetCenter(double e0, double e1, double e2)
    {
        center[0] = e0;
        center[1] = e1;
        center[2] = e2;
        Select(12);
    }

    public void SetRadius(double radius_)
    {
        radius = radius_;
        Select(13);
    }

    public void SetSphereInverse(boolean sphereInverse_)
    {
        sphereInverse = sphereInverse_;
        Select(14);
    }

    // Property getting methods
    public int      GetFuncType() { return funcType; }
    public boolean  GetPlane1Status() { return plane1Status; }
    public boolean  GetPlane2Status() { return plane2Status; }
    public boolean  GetPlane3Status() { return plane3Status; }
    public double[] GetPlane1Origin() { return plane1Origin; }
    public double[] GetPlane2Origin() { return plane2Origin; }
    public double[] GetPlane3Origin() { return plane3Origin; }
    public double[] GetPlane1Normal() { return plane1Normal; }
    public double[] GetPlane2Normal() { return plane2Normal; }
    public double[] GetPlane3Normal() { return plane3Normal; }
    public boolean  GetPlaneInverse() { return planeInverse; }
    public int      GetPlaneToolControlledClipPlane() { return planeToolControlledClipPlane; }
    public double[] GetCenter() { return center; }
    public double   GetRadius() { return radius; }
    public boolean  GetSphereInverse() { return sphereInverse; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(funcType);
        if(WriteSelect(1, buf))
            buf.WriteBool(plane1Status);
        if(WriteSelect(2, buf))
            buf.WriteBool(plane2Status);
        if(WriteSelect(3, buf))
            buf.WriteBool(plane3Status);
        if(WriteSelect(4, buf))
            buf.WriteDoubleArray(plane1Origin);
        if(WriteSelect(5, buf))
            buf.WriteDoubleArray(plane2Origin);
        if(WriteSelect(6, buf))
            buf.WriteDoubleArray(plane3Origin);
        if(WriteSelect(7, buf))
            buf.WriteDoubleArray(plane1Normal);
        if(WriteSelect(8, buf))
            buf.WriteDoubleArray(plane2Normal);
        if(WriteSelect(9, buf))
            buf.WriteDoubleArray(plane3Normal);
        if(WriteSelect(10, buf))
            buf.WriteBool(planeInverse);
        if(WriteSelect(11, buf))
            buf.WriteInt(planeToolControlledClipPlane);
        if(WriteSelect(12, buf))
            buf.WriteDoubleArray(center);
        if(WriteSelect(13, buf))
            buf.WriteDouble(radius);
        if(WriteSelect(14, buf))
            buf.WriteBool(sphereInverse);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetFuncType(buf.ReadInt());
                break;
            case 1:
                SetPlane1Status(buf.ReadBool());
                break;
            case 2:
                SetPlane2Status(buf.ReadBool());
                break;
            case 3:
                SetPlane3Status(buf.ReadBool());
                break;
            case 4:
                SetPlane1Origin(buf.ReadDoubleArray());
                break;
            case 5:
                SetPlane2Origin(buf.ReadDoubleArray());
                break;
            case 6:
                SetPlane3Origin(buf.ReadDoubleArray());
                break;
            case 7:
                SetPlane1Normal(buf.ReadDoubleArray());
                break;
            case 8:
                SetPlane2Normal(buf.ReadDoubleArray());
                break;
            case 9:
                SetPlane3Normal(buf.ReadDoubleArray());
                break;
            case 10:
                SetPlaneInverse(buf.ReadBool());
                break;
            case 11:
                SetPlaneToolControlledClipPlane(buf.ReadInt());
                break;
            case 12:
                SetCenter(buf.ReadDoubleArray());
                break;
            case 13:
                SetRadius(buf.ReadDouble());
                break;
            case 14:
                SetSphereInverse(buf.ReadBool());
                break;
            }
        }
    }


    // Attributes
    private int      funcType;
    private boolean  plane1Status;
    private boolean  plane2Status;
    private boolean  plane3Status;
    private double[] plane1Origin;
    private double[] plane2Origin;
    private double[] plane3Origin;
    private double[] plane1Normal;
    private double[] plane2Normal;
    private double[] plane3Normal;
    private boolean  planeInverse;
    private int      planeToolControlledClipPlane;
    private double[] center;
    private double   radius;
    private boolean  sphereInverse;
}

