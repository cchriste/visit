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

package llnl.visit.plots;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;
import llnl.visit.ColorAttribute;

// ****************************************************************************
// Class: PoincareAttributes
//
// Purpose:
//    Attributes for the Poincare plot
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class PoincareAttributes extends AttributeSubject implements Plugin
{
    // Enum values
    public final static int COLORINGMETHOD_SOLID = 0;
    public final static int COLORINGMETHOD_COLORBYSPEED = 1;
    public final static int COLORINGMETHOD_COLORBYVORTICITY = 2;
    public final static int COLORINGMETHOD_COLORBYLENGTH = 3;
    public final static int COLORINGMETHOD_COLORBYTIME = 4;
    public final static int COLORINGMETHOD_COLORBYSEEDPOINTID = 5;

    public final static int TERMINATIONTYPE_DISTANCE = 0;
    public final static int TERMINATIONTYPE_TIME = 1;


    public PoincareAttributes()
    {
        super(14);

        maxStepLength = 0.1;
        termination = 10;
        pointSource = new double[3];
        pointSource[0] = 0;
        pointSource[1] = 0;
        pointSource[2] = 0;
        planeOrigin = new double[3];
        planeOrigin[0] = 0;
        planeOrigin[1] = 0;
        planeOrigin[2] = 0;
        planeNormal = new double[3];
        planeNormal[0] = 0;
        planeNormal[1] = 0;
        planeNormal[2] = 1;
        planeUpAxis = new double[3];
        planeUpAxis[0] = 0;
        planeUpAxis[1] = 1;
        planeUpAxis[2] = 0;
        colorTableName = new String("Default");
        singleColor = new ColorAttribute(0, 0, 0);
        legendFlag = true;
        lightingFlag = true;
        relTol = 0.0001;
        absTol = 1e-05;
        terminationType = TERMINATIONTYPE_DISTANCE;
        integrationType = 0;
    }

    public PoincareAttributes(PoincareAttributes obj)
    {
        super(14);

        int i;

        maxStepLength = obj.maxStepLength;
        termination = obj.termination;
        pointSource = new double[3];
        pointSource[0] = obj.pointSource[0];
        pointSource[1] = obj.pointSource[1];
        pointSource[2] = obj.pointSource[2];

        planeOrigin = new double[3];
        planeOrigin[0] = obj.planeOrigin[0];
        planeOrigin[1] = obj.planeOrigin[1];
        planeOrigin[2] = obj.planeOrigin[2];

        planeNormal = new double[3];
        planeNormal[0] = obj.planeNormal[0];
        planeNormal[1] = obj.planeNormal[1];
        planeNormal[2] = obj.planeNormal[2];

        planeUpAxis = new double[3];
        planeUpAxis[0] = obj.planeUpAxis[0];
        planeUpAxis[1] = obj.planeUpAxis[1];
        planeUpAxis[2] = obj.planeUpAxis[2];

        colorTableName = new String(obj.colorTableName);
        singleColor = new ColorAttribute(obj.singleColor);
        legendFlag = obj.legendFlag;
        lightingFlag = obj.lightingFlag;
        relTol = obj.relTol;
        absTol = obj.absTol;
        terminationType = obj.terminationType;
        integrationType = obj.integrationType;

        SelectAll();
    }

    public boolean equals(PoincareAttributes obj)
    {
        int i;

        // Compare the pointSource arrays.
        boolean pointSource_equal = true;
        for(i = 0; i < 3 && pointSource_equal; ++i)
            pointSource_equal = (pointSource[i] == obj.pointSource[i]);

        // Compare the planeOrigin arrays.
        boolean planeOrigin_equal = true;
        for(i = 0; i < 3 && planeOrigin_equal; ++i)
            planeOrigin_equal = (planeOrigin[i] == obj.planeOrigin[i]);

        // Compare the planeNormal arrays.
        boolean planeNormal_equal = true;
        for(i = 0; i < 3 && planeNormal_equal; ++i)
            planeNormal_equal = (planeNormal[i] == obj.planeNormal[i]);

        // Compare the planeUpAxis arrays.
        boolean planeUpAxis_equal = true;
        for(i = 0; i < 3 && planeUpAxis_equal; ++i)
            planeUpAxis_equal = (planeUpAxis[i] == obj.planeUpAxis[i]);

        // Create the return value
        return ((maxStepLength == obj.maxStepLength) &&
                (termination == obj.termination) &&
                pointSource_equal &&
                planeOrigin_equal &&
                planeNormal_equal &&
                planeUpAxis_equal &&
                (colorTableName.equals(obj.colorTableName)) &&
                (singleColor == obj.singleColor) &&
                (legendFlag == obj.legendFlag) &&
                (lightingFlag == obj.lightingFlag) &&
                (relTol == obj.relTol) &&
                (absTol == obj.absTol) &&
                (terminationType == obj.terminationType) &&
                (integrationType == obj.integrationType));
    }

    public String GetName() { return "Poincare"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetMaxStepLength(double maxStepLength_)
    {
        maxStepLength = maxStepLength_;
        Select(0);
    }

    public void SetTermination(double termination_)
    {
        termination = termination_;
        Select(1);
    }

    public void SetPointSource(double[] pointSource_)
    {
        pointSource[0] = pointSource_[0];
        pointSource[1] = pointSource_[1];
        pointSource[2] = pointSource_[2];
        Select(2);
    }

    public void SetPointSource(double e0, double e1, double e2)
    {
        pointSource[0] = e0;
        pointSource[1] = e1;
        pointSource[2] = e2;
        Select(2);
    }

    public void SetPlaneOrigin(double[] planeOrigin_)
    {
        planeOrigin[0] = planeOrigin_[0];
        planeOrigin[1] = planeOrigin_[1];
        planeOrigin[2] = planeOrigin_[2];
        Select(3);
    }

    public void SetPlaneOrigin(double e0, double e1, double e2)
    {
        planeOrigin[0] = e0;
        planeOrigin[1] = e1;
        planeOrigin[2] = e2;
        Select(3);
    }

    public void SetPlaneNormal(double[] planeNormal_)
    {
        planeNormal[0] = planeNormal_[0];
        planeNormal[1] = planeNormal_[1];
        planeNormal[2] = planeNormal_[2];
        Select(4);
    }

    public void SetPlaneNormal(double e0, double e1, double e2)
    {
        planeNormal[0] = e0;
        planeNormal[1] = e1;
        planeNormal[2] = e2;
        Select(4);
    }

    public void SetPlaneUpAxis(double[] planeUpAxis_)
    {
        planeUpAxis[0] = planeUpAxis_[0];
        planeUpAxis[1] = planeUpAxis_[1];
        planeUpAxis[2] = planeUpAxis_[2];
        Select(5);
    }

    public void SetPlaneUpAxis(double e0, double e1, double e2)
    {
        planeUpAxis[0] = e0;
        planeUpAxis[1] = e1;
        planeUpAxis[2] = e2;
        Select(5);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(6);
    }

    public void SetSingleColor(ColorAttribute singleColor_)
    {
        singleColor = singleColor_;
        Select(7);
    }

    public void SetLegendFlag(boolean legendFlag_)
    {
        legendFlag = legendFlag_;
        Select(8);
    }

    public void SetLightingFlag(boolean lightingFlag_)
    {
        lightingFlag = lightingFlag_;
        Select(9);
    }

    public void SetRelTol(double relTol_)
    {
        relTol = relTol_;
        Select(10);
    }

    public void SetAbsTol(double absTol_)
    {
        absTol = absTol_;
        Select(11);
    }

    public void SetTerminationType(int terminationType_)
    {
        terminationType = terminationType_;
        Select(12);
    }

    public void SetIntegrationType(int integrationType_)
    {
        integrationType = integrationType_;
        Select(13);
    }

    // Property getting methods
    public double         GetMaxStepLength() { return maxStepLength; }
    public double         GetTermination() { return termination; }
    public double[]       GetPointSource() { return pointSource; }
    public double[]       GetPlaneOrigin() { return planeOrigin; }
    public double[]       GetPlaneNormal() { return planeNormal; }
    public double[]       GetPlaneUpAxis() { return planeUpAxis; }
    public String         GetColorTableName() { return colorTableName; }
    public ColorAttribute GetSingleColor() { return singleColor; }
    public boolean        GetLegendFlag() { return legendFlag; }
    public boolean        GetLightingFlag() { return lightingFlag; }
    public double         GetRelTol() { return relTol; }
    public double         GetAbsTol() { return absTol; }
    public int            GetTerminationType() { return terminationType; }
    public int            GetIntegrationType() { return integrationType; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDouble(maxStepLength);
        if(WriteSelect(1, buf))
            buf.WriteDouble(termination);
        if(WriteSelect(2, buf))
            buf.WriteDoubleArray(pointSource);
        if(WriteSelect(3, buf))
            buf.WriteDoubleArray(planeOrigin);
        if(WriteSelect(4, buf))
            buf.WriteDoubleArray(planeNormal);
        if(WriteSelect(5, buf))
            buf.WriteDoubleArray(planeUpAxis);
        if(WriteSelect(6, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(7, buf))
            singleColor.Write(buf);
        if(WriteSelect(8, buf))
            buf.WriteBool(legendFlag);
        if(WriteSelect(9, buf))
            buf.WriteBool(lightingFlag);
        if(WriteSelect(10, buf))
            buf.WriteDouble(relTol);
        if(WriteSelect(11, buf))
            buf.WriteDouble(absTol);
        if(WriteSelect(12, buf))
            buf.WriteInt(terminationType);
        if(WriteSelect(13, buf))
            buf.WriteInt(integrationType);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetMaxStepLength(buf.ReadDouble());
                break;
            case 1:
                SetTermination(buf.ReadDouble());
                break;
            case 2:
                SetPointSource(buf.ReadDoubleArray());
                break;
            case 3:
                SetPlaneOrigin(buf.ReadDoubleArray());
                break;
            case 4:
                SetPlaneNormal(buf.ReadDoubleArray());
                break;
            case 5:
                SetPlaneUpAxis(buf.ReadDoubleArray());
                break;
            case 6:
                SetColorTableName(buf.ReadString());
                break;
            case 7:
                singleColor.Read(buf);
                Select(7);
                break;
            case 8:
                SetLegendFlag(buf.ReadBool());
                break;
            case 9:
                SetLightingFlag(buf.ReadBool());
                break;
            case 10:
                SetRelTol(buf.ReadDouble());
                break;
            case 11:
                SetAbsTol(buf.ReadDouble());
                break;
            case 12:
                SetTerminationType(buf.ReadInt());
                break;
            case 13:
                SetIntegrationType(buf.ReadInt());
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleToString("maxStepLength", maxStepLength, indent) + "\n";
        str = str + doubleToString("termination", termination, indent) + "\n";
        str = str + doubleArrayToString("pointSource", pointSource, indent) + "\n";
        str = str + doubleArrayToString("planeOrigin", planeOrigin, indent) + "\n";
        str = str + doubleArrayToString("planeNormal", planeNormal, indent) + "\n";
        str = str + doubleArrayToString("planeUpAxis", planeUpAxis, indent) + "\n";
        str = str + stringToString("colorTableName", colorTableName, indent) + "\n";
        str = str + indent + "singleColor = {" + singleColor.Red() + ", " + singleColor.Green() + ", " + singleColor.Blue() + ", " + singleColor.Alpha() + "}\n";
        str = str + boolToString("legendFlag", legendFlag, indent) + "\n";
        str = str + boolToString("lightingFlag", lightingFlag, indent) + "\n";
        str = str + doubleToString("relTol", relTol, indent) + "\n";
        str = str + doubleToString("absTol", absTol, indent) + "\n";
        str = str + indent + "terminationType = ";
        if(terminationType == TERMINATIONTYPE_DISTANCE)
            str = str + "TERMINATIONTYPE_DISTANCE";
        if(terminationType == TERMINATIONTYPE_TIME)
            str = str + "TERMINATIONTYPE_TIME";
        str = str + "\n";
        str = str + intToString("integrationType", integrationType, indent) + "\n";
        return str;
    }


    // Attributes
    private double         maxStepLength;
    private double         termination;
    private double[]       pointSource;
    private double[]       planeOrigin;
    private double[]       planeNormal;
    private double[]       planeUpAxis;
    private String         colorTableName;
    private ColorAttribute singleColor;
    private boolean        legendFlag;
    private boolean        lightingFlag;
    private double         relTol;
    private double         absTol;
    private int            terminationType;
    private int            integrationType;
}

