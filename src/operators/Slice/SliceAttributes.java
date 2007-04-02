// ****************************************************************************
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
// ****************************************************************************

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: SliceAttributes
//
// Purpose:
//    This class contains attributes for the arbitrary slice.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Tue Jan 25 16:57:14 PST 2005
//
// Modifications:
//   
// ****************************************************************************

public class SliceAttributes extends AttributeSubject implements Plugin
{
    // Enum values
    public final static int AXISTYPE_XAXIS = 0;
    public final static int AXISTYPE_YAXIS = 1;
    public final static int AXISTYPE_ZAXIS = 2;
    public final static int AXISTYPE_ARBITRARY = 3;

    public final static int ORIGINTYPE_POINT = 0;
    public final static int ORIGINTYPE_INTERCEPT = 1;
    public final static int ORIGINTYPE_PERCENT = 2;
    public final static int ORIGINTYPE_ZONE = 3;
    public final static int ORIGINTYPE_NODE = 4;


    public SliceAttributes()
    {
        super(15);

        originType = ORIGINTYPE_INTERCEPT;
        originPoint = new double[3];
        originPoint[0] = 0;
        originPoint[1] = 0;
        originPoint[2] = 0;
        originIntercept = 0;
        originPercent = 0;
        originZone = 0;
        originNode = 0;
        normal = new double[3];
        normal[0] = 0;
        normal[1] = -1;
        normal[2] = 0;
        axisType = AXISTYPE_YAXIS;
        upAxis = new double[3];
        upAxis[0] = 0;
        upAxis[1] = 0;
        upAxis[2] = 1;
        project2d = true;
        interactive = true;
        flip = false;
        originZoneDomain = 0;
        originNodeDomain = 0;
        meshName = new String("default");
    }

    public SliceAttributes(SliceAttributes obj)
    {
        super(15);

        int i;

        originType = obj.originType;
        originPoint = new double[3];
        originPoint[0] = obj.originPoint[0];
        originPoint[1] = obj.originPoint[1];
        originPoint[2] = obj.originPoint[2];

        originIntercept = obj.originIntercept;
        originPercent = obj.originPercent;
        originZone = obj.originZone;
        originNode = obj.originNode;
        normal = new double[3];
        normal[0] = obj.normal[0];
        normal[1] = obj.normal[1];
        normal[2] = obj.normal[2];

        axisType = obj.axisType;
        upAxis = new double[3];
        upAxis[0] = obj.upAxis[0];
        upAxis[1] = obj.upAxis[1];
        upAxis[2] = obj.upAxis[2];

        project2d = obj.project2d;
        interactive = obj.interactive;
        flip = obj.flip;
        originZoneDomain = obj.originZoneDomain;
        originNodeDomain = obj.originNodeDomain;
        meshName = new String(obj.meshName);

        SelectAll();
    }

    public boolean equals(SliceAttributes obj)
    {
        int i;

        // Compare the originPoint arrays.
        boolean originPoint_equal = true;
        for(i = 0; i < 3 && originPoint_equal; ++i)
            originPoint_equal = (originPoint[i] == obj.originPoint[i]);

        // Compare the normal arrays.
        boolean normal_equal = true;
        for(i = 0; i < 3 && normal_equal; ++i)
            normal_equal = (normal[i] == obj.normal[i]);

        // Compare the upAxis arrays.
        boolean upAxis_equal = true;
        for(i = 0; i < 3 && upAxis_equal; ++i)
            upAxis_equal = (upAxis[i] == obj.upAxis[i]);

        // Create the return value
        return ((originType == obj.originType) &&
                originPoint_equal &&
                (originIntercept == obj.originIntercept) &&
                (originPercent == obj.originPercent) &&
                (originZone == obj.originZone) &&
                (originNode == obj.originNode) &&
                normal_equal &&
                true /* can ignore axisType */ &&
                upAxis_equal &&
                (project2d == obj.project2d) &&
                (interactive == obj.interactive) &&
                true /* can ignore flip */ &&
                (originZoneDomain == obj.originZoneDomain) &&
                (originNodeDomain == obj.originNodeDomain) &&
                (meshName == obj.meshName));
    }

    public String GetName() { return "Slice"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetOriginType(int originType_)
    {
        originType = originType_;
        Select(0);
    }

    public void SetOriginPoint(double[] originPoint_)
    {
        originPoint[0] = originPoint_[0];
        originPoint[1] = originPoint_[1];
        originPoint[2] = originPoint_[2];
        Select(1);
    }

    public void SetOriginPoint(double e0, double e1, double e2)
    {
        originPoint[0] = e0;
        originPoint[1] = e1;
        originPoint[2] = e2;
        Select(1);
    }

    public void SetOriginIntercept(double originIntercept_)
    {
        originIntercept = originIntercept_;
        Select(2);
    }

    public void SetOriginPercent(double originPercent_)
    {
        originPercent = originPercent_;
        Select(3);
    }

    public void SetOriginZone(int originZone_)
    {
        originZone = originZone_;
        Select(4);
    }

    public void SetOriginNode(int originNode_)
    {
        originNode = originNode_;
        Select(5);
    }

    public void SetNormal(double[] normal_)
    {
        normal[0] = normal_[0];
        normal[1] = normal_[1];
        normal[2] = normal_[2];
        Select(6);
    }

    public void SetNormal(double e0, double e1, double e2)
    {
        normal[0] = e0;
        normal[1] = e1;
        normal[2] = e2;
        Select(6);
    }

    public void SetAxisType(int axisType_)
    {
        axisType = axisType_;
        Select(7);
    }

    public void SetUpAxis(double[] upAxis_)
    {
        upAxis[0] = upAxis_[0];
        upAxis[1] = upAxis_[1];
        upAxis[2] = upAxis_[2];
        Select(8);
    }

    public void SetUpAxis(double e0, double e1, double e2)
    {
        upAxis[0] = e0;
        upAxis[1] = e1;
        upAxis[2] = e2;
        Select(8);
    }

    public void SetProject2d(boolean project2d_)
    {
        project2d = project2d_;
        Select(9);
    }

    public void SetInteractive(boolean interactive_)
    {
        interactive = interactive_;
        Select(10);
    }

    public void SetFlip(boolean flip_)
    {
        flip = flip_;
        Select(11);
    }

    public void SetOriginZoneDomain(int originZoneDomain_)
    {
        originZoneDomain = originZoneDomain_;
        Select(12);
    }

    public void SetOriginNodeDomain(int originNodeDomain_)
    {
        originNodeDomain = originNodeDomain_;
        Select(13);
    }

    public void SetMeshName(String meshName_)
    {
        meshName = meshName_;
        Select(14);
    }

    // Property getting methods
    public int      GetOriginType() { return originType; }
    public double[] GetOriginPoint() { return originPoint; }
    public double   GetOriginIntercept() { return originIntercept; }
    public double   GetOriginPercent() { return originPercent; }
    public int      GetOriginZone() { return originZone; }
    public int      GetOriginNode() { return originNode; }
    public double[] GetNormal() { return normal; }
    public int      GetAxisType() { return axisType; }
    public double[] GetUpAxis() { return upAxis; }
    public boolean  GetProject2d() { return project2d; }
    public boolean  GetInteractive() { return interactive; }
    public boolean  GetFlip() { return flip; }
    public int      GetOriginZoneDomain() { return originZoneDomain; }
    public int      GetOriginNodeDomain() { return originNodeDomain; }
    public String   GetMeshName() { return meshName; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(originType);
        if(WriteSelect(1, buf))
            buf.WriteDoubleArray(originPoint);
        if(WriteSelect(2, buf))
            buf.WriteDouble(originIntercept);
        if(WriteSelect(3, buf))
            buf.WriteDouble(originPercent);
        if(WriteSelect(4, buf))
            buf.WriteInt(originZone);
        if(WriteSelect(5, buf))
            buf.WriteInt(originNode);
        if(WriteSelect(6, buf))
            buf.WriteDoubleArray(normal);
        if(WriteSelect(7, buf))
            buf.WriteInt(axisType);
        if(WriteSelect(8, buf))
            buf.WriteDoubleArray(upAxis);
        if(WriteSelect(9, buf))
            buf.WriteBool(project2d);
        if(WriteSelect(10, buf))
            buf.WriteBool(interactive);
        if(WriteSelect(11, buf))
            buf.WriteBool(flip);
        if(WriteSelect(12, buf))
            buf.WriteInt(originZoneDomain);
        if(WriteSelect(13, buf))
            buf.WriteInt(originNodeDomain);
        if(WriteSelect(14, buf))
            buf.WriteString(meshName);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetOriginType(buf.ReadInt());
                break;
            case 1:
                SetOriginPoint(buf.ReadDoubleArray());
                break;
            case 2:
                SetOriginIntercept(buf.ReadDouble());
                break;
            case 3:
                SetOriginPercent(buf.ReadDouble());
                break;
            case 4:
                SetOriginZone(buf.ReadInt());
                break;
            case 5:
                SetOriginNode(buf.ReadInt());
                break;
            case 6:
                SetNormal(buf.ReadDoubleArray());
                break;
            case 7:
                SetAxisType(buf.ReadInt());
                break;
            case 8:
                SetUpAxis(buf.ReadDoubleArray());
                break;
            case 9:
                SetProject2d(buf.ReadBool());
                break;
            case 10:
                SetInteractive(buf.ReadBool());
                break;
            case 11:
                SetFlip(buf.ReadBool());
                break;
            case 12:
                SetOriginZoneDomain(buf.ReadInt());
                break;
            case 13:
                SetOriginNodeDomain(buf.ReadInt());
                break;
            case 14:
                SetMeshName(buf.ReadString());
                break;
            }
        }
    }


    // Attributes
    private int      originType;
    private double[] originPoint;
    private double   originIntercept;
    private double   originPercent;
    private int      originZone;
    private int      originNode;
    private double[] normal;
    private int      axisType;
    private double[] upAxis;
    private boolean  project2d;
    private boolean  interactive;
    private boolean  flip;
    private int      originZoneDomain;
    private int      originNodeDomain;
    private String   meshName;
}

