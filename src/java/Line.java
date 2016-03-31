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

package llnl.visit;


// ****************************************************************************
// Class: Line
//
// Purpose:
//    Attributes for a line
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class Line extends AttributeSubject
{
    private static int Line_numAdditionalAtts = 11;

    public Line()
    {
        super(Line_numAdditionalAtts);

        point1 = new double[3];
        point1[0] = 0;
        point1[1] = 0;
        point1[2] = 0;
        point2 = new double[3];
        point2[0] = 1;
        point2[1] = 1;
        point2[2] = 1;
        designator = new String("");
        color = new ColorAttribute(0, 0, 0);
        lineStyle = 0;
        lineWidth = 0;
        interactive = false;
        reflineLabels = false;
        varName = new String("default");
        samplingOn = false;
        numSamples = 50;
    }

    public Line(int nMoreFields)
    {
        super(Line_numAdditionalAtts + nMoreFields);

        point1 = new double[3];
        point1[0] = 0;
        point1[1] = 0;
        point1[2] = 0;
        point2 = new double[3];
        point2[0] = 1;
        point2[1] = 1;
        point2[2] = 1;
        designator = new String("");
        color = new ColorAttribute(0, 0, 0);
        lineStyle = 0;
        lineWidth = 0;
        interactive = false;
        reflineLabels = false;
        varName = new String("default");
        samplingOn = false;
        numSamples = 50;
    }

    public Line(Line obj)
    {
        super(Line_numAdditionalAtts);

        int i;

        point1 = new double[3];
        point1[0] = obj.point1[0];
        point1[1] = obj.point1[1];
        point1[2] = obj.point1[2];

        point2 = new double[3];
        point2[0] = obj.point2[0];
        point2[1] = obj.point2[1];
        point2[2] = obj.point2[2];

        designator = new String(obj.designator);
        color = new ColorAttribute(obj.color);
        lineStyle = obj.lineStyle;
        lineWidth = obj.lineWidth;
        interactive = obj.interactive;
        reflineLabels = obj.reflineLabels;
        varName = new String(obj.varName);
        samplingOn = obj.samplingOn;
        numSamples = obj.numSamples;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return Line_numAdditionalAtts;
    }

    public boolean equals(Line obj)
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
                (designator.equals(obj.designator)) &&
                (color == obj.color) &&
                (lineStyle == obj.lineStyle) &&
                (lineWidth == obj.lineWidth) &&
                (interactive == obj.interactive) &&
                (reflineLabels == obj.reflineLabels) &&
                (varName.equals(obj.varName)) &&
                (samplingOn == obj.samplingOn) &&
                (numSamples == obj.numSamples));
    }

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

    public void SetDesignator(String designator_)
    {
        designator = designator_;
        Select(2);
    }

    public void SetColor(ColorAttribute color_)
    {
        color = color_;
        Select(3);
    }

    public void SetLineStyle(int lineStyle_)
    {
        lineStyle = lineStyle_;
        Select(4);
    }

    public void SetLineWidth(int lineWidth_)
    {
        lineWidth = lineWidth_;
        Select(5);
    }

    public void SetInteractive(boolean interactive_)
    {
        interactive = interactive_;
        Select(6);
    }

    public void SetReflineLabels(boolean reflineLabels_)
    {
        reflineLabels = reflineLabels_;
        Select(7);
    }

    public void SetVarName(String varName_)
    {
        varName = varName_;
        Select(8);
    }

    public void SetSamplingOn(boolean samplingOn_)
    {
        samplingOn = samplingOn_;
        Select(9);
    }

    public void SetNumSamples(int numSamples_)
    {
        numSamples = numSamples_;
        Select(10);
    }

    // Property getting methods
    public double[]       GetPoint1() { return point1; }
    public double[]       GetPoint2() { return point2; }
    public String         GetDesignator() { return designator; }
    public ColorAttribute GetColor() { return color; }
    public int            GetLineStyle() { return lineStyle; }
    public int            GetLineWidth() { return lineWidth; }
    public boolean        GetInteractive() { return interactive; }
    public boolean        GetReflineLabels() { return reflineLabels; }
    public String         GetVarName() { return varName; }
    public boolean        GetSamplingOn() { return samplingOn; }
    public int            GetNumSamples() { return numSamples; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(point1);
        if(WriteSelect(1, buf))
            buf.WriteDoubleArray(point2);
        if(WriteSelect(2, buf))
            buf.WriteString(designator);
        if(WriteSelect(3, buf))
            color.Write(buf);
        if(WriteSelect(4, buf))
            buf.WriteInt(lineStyle);
        if(WriteSelect(5, buf))
            buf.WriteInt(lineWidth);
        if(WriteSelect(6, buf))
            buf.WriteBool(interactive);
        if(WriteSelect(7, buf))
            buf.WriteBool(reflineLabels);
        if(WriteSelect(8, buf))
            buf.WriteString(varName);
        if(WriteSelect(9, buf))
            buf.WriteBool(samplingOn);
        if(WriteSelect(10, buf))
            buf.WriteInt(numSamples);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetPoint1(buf.ReadDoubleArray());
            break;
        case 1:
            SetPoint2(buf.ReadDoubleArray());
            break;
        case 2:
            SetDesignator(buf.ReadString());
            break;
        case 3:
            color.Read(buf);
            Select(3);
            break;
        case 4:
            SetLineStyle(buf.ReadInt());
            break;
        case 5:
            SetLineWidth(buf.ReadInt());
            break;
        case 6:
            SetInteractive(buf.ReadBool());
            break;
        case 7:
            SetReflineLabels(buf.ReadBool());
            break;
        case 8:
            SetVarName(buf.ReadString());
            break;
        case 9:
            SetSamplingOn(buf.ReadBool());
            break;
        case 10:
            SetNumSamples(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + doubleArrayToString("point1", point1, indent) + "\n";
        str = str + doubleArrayToString("point2", point2, indent) + "\n";
        str = str + stringToString("designator", designator, indent) + "\n";
        str = str + indent + "color = {" + color.Red() + ", " + color.Green() + ", " + color.Blue() + ", " + color.Alpha() + "}\n";
        str = str + intToString("lineStyle", lineStyle, indent) + "\n";
        str = str + intToString("lineWidth", lineWidth, indent) + "\n";
        str = str + boolToString("interactive", interactive, indent) + "\n";
        str = str + boolToString("reflineLabels", reflineLabels, indent) + "\n";
        str = str + stringToString("varName", varName, indent) + "\n";
        str = str + boolToString("samplingOn", samplingOn, indent) + "\n";
        str = str + intToString("numSamples", numSamples, indent) + "\n";
        return str;
    }


    // Attributes
    private double[]       point1;
    private double[]       point2;
    private String         designator;
    private ColorAttribute color;
    private int            lineStyle;
    private int            lineWidth;
    private boolean        interactive;
    private boolean        reflineLabels;
    private String         varName;
    private boolean        samplingOn;
    private int            numSamples;
}

