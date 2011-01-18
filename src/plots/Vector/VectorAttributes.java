// ***************************************************************************
//
// Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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

package llnl.visit.plots;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;
import llnl.visit.ColorAttribute;

// ****************************************************************************
// Class: VectorAttributes
//
// Purpose:
//    Attributes for the vector plot
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class VectorAttributes extends AttributeSubject implements Plugin
{
    private static int VectorAttributes_numAdditionalAtts = 27;

    // Enum values
    public final static int QUALITY_FAST = 0;
    public final static int QUALITY_HIGH = 1;

    public final static int ORIGINTYPE_HEAD = 0;
    public final static int ORIGINTYPE_MIDDLE = 1;
    public final static int ORIGINTYPE_TAIL = 2;

    public final static int LIMITSMODE_ORIGINALDATA = 0;
    public final static int LIMITSMODE_CURRENTPLOT = 1;

    public final static int GLYPHTYPE_ARROW = 0;
    public final static int GLYPHTYPE_ELLIPSOID = 1;

    public final static int GLYPHLOCATION_ADAPTSTOMESHRESOLUTION = 0;
    public final static int GLYPHLOCATION_UNIFORMINSPACE = 1;


    public VectorAttributes()
    {
        super(VectorAttributes_numAdditionalAtts);

        glyphLocation = GLYPHLOCATION_ADAPTSTOMESHRESOLUTION;
        useStride = false;
        stride = 1;
        nVectors = 400;
        lineStyle = 0;
        lineWidth = 0;
        scale = 0.25;
        scaleByMagnitude = true;
        autoScale = true;
        headSize = 0.25;
        headOn = true;
        colorByMag = true;
        useLegend = true;
        vectorColor = new ColorAttribute(0, 0, 0);
        colorTableName = new String("Default");
        invertColorTable = false;
        vectorOrigin = ORIGINTYPE_TAIL;
        minFlag = false;
        maxFlag = false;
        limitsMode = LIMITSMODE_ORIGINALDATA;
        min = 0;
        max = 1;
        lineStem = true;
        geometryQuality = QUALITY_FAST;
        stemWidth = 0.08;
        origOnly = true;
        glyphType = GLYPHTYPE_ARROW;
    }

    public VectorAttributes(int nMoreFields)
    {
        super(VectorAttributes_numAdditionalAtts + nMoreFields);

        glyphLocation = GLYPHLOCATION_ADAPTSTOMESHRESOLUTION;
        useStride = false;
        stride = 1;
        nVectors = 400;
        lineStyle = 0;
        lineWidth = 0;
        scale = 0.25;
        scaleByMagnitude = true;
        autoScale = true;
        headSize = 0.25;
        headOn = true;
        colorByMag = true;
        useLegend = true;
        vectorColor = new ColorAttribute(0, 0, 0);
        colorTableName = new String("Default");
        invertColorTable = false;
        vectorOrigin = ORIGINTYPE_TAIL;
        minFlag = false;
        maxFlag = false;
        limitsMode = LIMITSMODE_ORIGINALDATA;
        min = 0;
        max = 1;
        lineStem = true;
        geometryQuality = QUALITY_FAST;
        stemWidth = 0.08;
        origOnly = true;
        glyphType = GLYPHTYPE_ARROW;
    }

    public VectorAttributes(VectorAttributes obj)
    {
        super(VectorAttributes_numAdditionalAtts);

        glyphLocation = obj.glyphLocation;
        useStride = obj.useStride;
        stride = obj.stride;
        nVectors = obj.nVectors;
        lineStyle = obj.lineStyle;
        lineWidth = obj.lineWidth;
        scale = obj.scale;
        scaleByMagnitude = obj.scaleByMagnitude;
        autoScale = obj.autoScale;
        headSize = obj.headSize;
        headOn = obj.headOn;
        colorByMag = obj.colorByMag;
        useLegend = obj.useLegend;
        vectorColor = new ColorAttribute(obj.vectorColor);
        colorTableName = new String(obj.colorTableName);
        invertColorTable = obj.invertColorTable;
        vectorOrigin = obj.vectorOrigin;
        minFlag = obj.minFlag;
        maxFlag = obj.maxFlag;
        limitsMode = obj.limitsMode;
        min = obj.min;
        max = obj.max;
        lineStem = obj.lineStem;
        geometryQuality = obj.geometryQuality;
        stemWidth = obj.stemWidth;
        origOnly = obj.origOnly;
        glyphType = obj.glyphType;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return VectorAttributes_numAdditionalAtts;
    }

    public boolean equals(VectorAttributes obj)
    {
        // Create the return value
        return ((glyphLocation == obj.glyphLocation) &&
                (useStride == obj.useStride) &&
                (stride == obj.stride) &&
                (nVectors == obj.nVectors) &&
                (lineStyle == obj.lineStyle) &&
                (lineWidth == obj.lineWidth) &&
                (scale == obj.scale) &&
                (scaleByMagnitude == obj.scaleByMagnitude) &&
                (autoScale == obj.autoScale) &&
                (headSize == obj.headSize) &&
                (headOn == obj.headOn) &&
                (colorByMag == obj.colorByMag) &&
                (useLegend == obj.useLegend) &&
                (vectorColor == obj.vectorColor) &&
                (colorTableName.equals(obj.colorTableName)) &&
                (invertColorTable == obj.invertColorTable) &&
                (vectorOrigin == obj.vectorOrigin) &&
                (minFlag == obj.minFlag) &&
                (maxFlag == obj.maxFlag) &&
                (limitsMode == obj.limitsMode) &&
                (min == obj.min) &&
                (max == obj.max) &&
                (lineStem == obj.lineStem) &&
                (geometryQuality == obj.geometryQuality) &&
                (stemWidth == obj.stemWidth) &&
                (origOnly == obj.origOnly) &&
                (glyphType == obj.glyphType));
    }

    public String GetName() { return "Vector"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetGlyphLocation(int glyphLocation_)
    {
        glyphLocation = glyphLocation_;
        Select(0);
    }

    public void SetUseStride(boolean useStride_)
    {
        useStride = useStride_;
        Select(1);
    }

    public void SetStride(int stride_)
    {
        stride = stride_;
        Select(2);
    }

    public void SetNVectors(int nVectors_)
    {
        nVectors = nVectors_;
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

    public void SetScale(double scale_)
    {
        scale = scale_;
        Select(6);
    }

    public void SetScaleByMagnitude(boolean scaleByMagnitude_)
    {
        scaleByMagnitude = scaleByMagnitude_;
        Select(7);
    }

    public void SetAutoScale(boolean autoScale_)
    {
        autoScale = autoScale_;
        Select(8);
    }

    public void SetHeadSize(double headSize_)
    {
        headSize = headSize_;
        Select(9);
    }

    public void SetHeadOn(boolean headOn_)
    {
        headOn = headOn_;
        Select(10);
    }

    public void SetColorByMag(boolean colorByMag_)
    {
        colorByMag = colorByMag_;
        Select(11);
    }

    public void SetUseLegend(boolean useLegend_)
    {
        useLegend = useLegend_;
        Select(12);
    }

    public void SetVectorColor(ColorAttribute vectorColor_)
    {
        vectorColor = vectorColor_;
        Select(13);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(14);
    }

    public void SetInvertColorTable(boolean invertColorTable_)
    {
        invertColorTable = invertColorTable_;
        Select(15);
    }

    public void SetVectorOrigin(int vectorOrigin_)
    {
        vectorOrigin = vectorOrigin_;
        Select(16);
    }

    public void SetMinFlag(boolean minFlag_)
    {
        minFlag = minFlag_;
        Select(17);
    }

    public void SetMaxFlag(boolean maxFlag_)
    {
        maxFlag = maxFlag_;
        Select(18);
    }

    public void SetLimitsMode(int limitsMode_)
    {
        limitsMode = limitsMode_;
        Select(19);
    }

    public void SetMin(double min_)
    {
        min = min_;
        Select(20);
    }

    public void SetMax(double max_)
    {
        max = max_;
        Select(21);
    }

    public void SetLineStem(boolean lineStem_)
    {
        lineStem = lineStem_;
        Select(22);
    }

    public void SetGeometryQuality(int geometryQuality_)
    {
        geometryQuality = geometryQuality_;
        Select(23);
    }

    public void SetStemWidth(double stemWidth_)
    {
        stemWidth = stemWidth_;
        Select(24);
    }

    public void SetOrigOnly(boolean origOnly_)
    {
        origOnly = origOnly_;
        Select(25);
    }

    public void SetGlyphType(int glyphType_)
    {
        glyphType = glyphType_;
        Select(26);
    }

    // Property getting methods
    public int            GetGlyphLocation() { return glyphLocation; }
    public boolean        GetUseStride() { return useStride; }
    public int            GetStride() { return stride; }
    public int            GetNVectors() { return nVectors; }
    public int            GetLineStyle() { return lineStyle; }
    public int            GetLineWidth() { return lineWidth; }
    public double         GetScale() { return scale; }
    public boolean        GetScaleByMagnitude() { return scaleByMagnitude; }
    public boolean        GetAutoScale() { return autoScale; }
    public double         GetHeadSize() { return headSize; }
    public boolean        GetHeadOn() { return headOn; }
    public boolean        GetColorByMag() { return colorByMag; }
    public boolean        GetUseLegend() { return useLegend; }
    public ColorAttribute GetVectorColor() { return vectorColor; }
    public String         GetColorTableName() { return colorTableName; }
    public boolean        GetInvertColorTable() { return invertColorTable; }
    public int            GetVectorOrigin() { return vectorOrigin; }
    public boolean        GetMinFlag() { return minFlag; }
    public boolean        GetMaxFlag() { return maxFlag; }
    public int            GetLimitsMode() { return limitsMode; }
    public double         GetMin() { return min; }
    public double         GetMax() { return max; }
    public boolean        GetLineStem() { return lineStem; }
    public int            GetGeometryQuality() { return geometryQuality; }
    public double         GetStemWidth() { return stemWidth; }
    public boolean        GetOrigOnly() { return origOnly; }
    public int            GetGlyphType() { return glyphType; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(glyphLocation);
        if(WriteSelect(1, buf))
            buf.WriteBool(useStride);
        if(WriteSelect(2, buf))
            buf.WriteInt(stride);
        if(WriteSelect(3, buf))
            buf.WriteInt(nVectors);
        if(WriteSelect(4, buf))
            buf.WriteInt(lineStyle);
        if(WriteSelect(5, buf))
            buf.WriteInt(lineWidth);
        if(WriteSelect(6, buf))
            buf.WriteDouble(scale);
        if(WriteSelect(7, buf))
            buf.WriteBool(scaleByMagnitude);
        if(WriteSelect(8, buf))
            buf.WriteBool(autoScale);
        if(WriteSelect(9, buf))
            buf.WriteDouble(headSize);
        if(WriteSelect(10, buf))
            buf.WriteBool(headOn);
        if(WriteSelect(11, buf))
            buf.WriteBool(colorByMag);
        if(WriteSelect(12, buf))
            buf.WriteBool(useLegend);
        if(WriteSelect(13, buf))
            vectorColor.Write(buf);
        if(WriteSelect(14, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(15, buf))
            buf.WriteBool(invertColorTable);
        if(WriteSelect(16, buf))
            buf.WriteInt(vectorOrigin);
        if(WriteSelect(17, buf))
            buf.WriteBool(minFlag);
        if(WriteSelect(18, buf))
            buf.WriteBool(maxFlag);
        if(WriteSelect(19, buf))
            buf.WriteInt(limitsMode);
        if(WriteSelect(20, buf))
            buf.WriteDouble(min);
        if(WriteSelect(21, buf))
            buf.WriteDouble(max);
        if(WriteSelect(22, buf))
            buf.WriteBool(lineStem);
        if(WriteSelect(23, buf))
            buf.WriteInt(geometryQuality);
        if(WriteSelect(24, buf))
            buf.WriteDouble(stemWidth);
        if(WriteSelect(25, buf))
            buf.WriteBool(origOnly);
        if(WriteSelect(26, buf))
            buf.WriteInt(glyphType);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetGlyphLocation(buf.ReadInt());
            break;
        case 1:
            SetUseStride(buf.ReadBool());
            break;
        case 2:
            SetStride(buf.ReadInt());
            break;
        case 3:
            SetNVectors(buf.ReadInt());
            break;
        case 4:
            SetLineStyle(buf.ReadInt());
            break;
        case 5:
            SetLineWidth(buf.ReadInt());
            break;
        case 6:
            SetScale(buf.ReadDouble());
            break;
        case 7:
            SetScaleByMagnitude(buf.ReadBool());
            break;
        case 8:
            SetAutoScale(buf.ReadBool());
            break;
        case 9:
            SetHeadSize(buf.ReadDouble());
            break;
        case 10:
            SetHeadOn(buf.ReadBool());
            break;
        case 11:
            SetColorByMag(buf.ReadBool());
            break;
        case 12:
            SetUseLegend(buf.ReadBool());
            break;
        case 13:
            vectorColor.Read(buf);
            Select(13);
            break;
        case 14:
            SetColorTableName(buf.ReadString());
            break;
        case 15:
            SetInvertColorTable(buf.ReadBool());
            break;
        case 16:
            SetVectorOrigin(buf.ReadInt());
            break;
        case 17:
            SetMinFlag(buf.ReadBool());
            break;
        case 18:
            SetMaxFlag(buf.ReadBool());
            break;
        case 19:
            SetLimitsMode(buf.ReadInt());
            break;
        case 20:
            SetMin(buf.ReadDouble());
            break;
        case 21:
            SetMax(buf.ReadDouble());
            break;
        case 22:
            SetLineStem(buf.ReadBool());
            break;
        case 23:
            SetGeometryQuality(buf.ReadInt());
            break;
        case 24:
            SetStemWidth(buf.ReadDouble());
            break;
        case 25:
            SetOrigOnly(buf.ReadBool());
            break;
        case 26:
            SetGlyphType(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "glyphLocation = ";
        if(glyphLocation == GLYPHLOCATION_ADAPTSTOMESHRESOLUTION)
            str = str + "GLYPHLOCATION_ADAPTSTOMESHRESOLUTION";
        if(glyphLocation == GLYPHLOCATION_UNIFORMINSPACE)
            str = str + "GLYPHLOCATION_UNIFORMINSPACE";
        str = str + "\n";
        str = str + boolToString("useStride", useStride, indent) + "\n";
        str = str + intToString("stride", stride, indent) + "\n";
        str = str + intToString("nVectors", nVectors, indent) + "\n";
        str = str + intToString("lineStyle", lineStyle, indent) + "\n";
        str = str + intToString("lineWidth", lineWidth, indent) + "\n";
        str = str + doubleToString("scale", scale, indent) + "\n";
        str = str + boolToString("scaleByMagnitude", scaleByMagnitude, indent) + "\n";
        str = str + boolToString("autoScale", autoScale, indent) + "\n";
        str = str + doubleToString("headSize", headSize, indent) + "\n";
        str = str + boolToString("headOn", headOn, indent) + "\n";
        str = str + boolToString("colorByMag", colorByMag, indent) + "\n";
        str = str + boolToString("useLegend", useLegend, indent) + "\n";
        str = str + indent + "vectorColor = {" + vectorColor.Red() + ", " + vectorColor.Green() + ", " + vectorColor.Blue() + ", " + vectorColor.Alpha() + "}\n";
        str = str + stringToString("colorTableName", colorTableName, indent) + "\n";
        str = str + boolToString("invertColorTable", invertColorTable, indent) + "\n";
        str = str + indent + "vectorOrigin = ";
        if(vectorOrigin == ORIGINTYPE_HEAD)
            str = str + "ORIGINTYPE_HEAD";
        if(vectorOrigin == ORIGINTYPE_MIDDLE)
            str = str + "ORIGINTYPE_MIDDLE";
        if(vectorOrigin == ORIGINTYPE_TAIL)
            str = str + "ORIGINTYPE_TAIL";
        str = str + "\n";
        str = str + boolToString("minFlag", minFlag, indent) + "\n";
        str = str + boolToString("maxFlag", maxFlag, indent) + "\n";
        str = str + indent + "limitsMode = ";
        if(limitsMode == LIMITSMODE_ORIGINALDATA)
            str = str + "LIMITSMODE_ORIGINALDATA";
        if(limitsMode == LIMITSMODE_CURRENTPLOT)
            str = str + "LIMITSMODE_CURRENTPLOT";
        str = str + "\n";
        str = str + doubleToString("min", min, indent) + "\n";
        str = str + doubleToString("max", max, indent) + "\n";
        str = str + boolToString("lineStem", lineStem, indent) + "\n";
        str = str + indent + "geometryQuality = ";
        if(geometryQuality == QUALITY_FAST)
            str = str + "QUALITY_FAST";
        if(geometryQuality == QUALITY_HIGH)
            str = str + "QUALITY_HIGH";
        str = str + "\n";
        str = str + doubleToString("stemWidth", stemWidth, indent) + "\n";
        str = str + boolToString("origOnly", origOnly, indent) + "\n";
        str = str + indent + "glyphType = ";
        if(glyphType == GLYPHTYPE_ARROW)
            str = str + "GLYPHTYPE_ARROW";
        if(glyphType == GLYPHTYPE_ELLIPSOID)
            str = str + "GLYPHTYPE_ELLIPSOID";
        str = str + "\n";
        return str;
    }


    // Attributes
    private int            glyphLocation;
    private boolean        useStride;
    private int            stride;
    private int            nVectors;
    private int            lineStyle;
    private int            lineWidth;
    private double         scale;
    private boolean        scaleByMagnitude;
    private boolean        autoScale;
    private double         headSize;
    private boolean        headOn;
    private boolean        colorByMag;
    private boolean        useLegend;
    private ColorAttribute vectorColor;
    private String         colorTableName;
    private boolean        invertColorTable;
    private int            vectorOrigin;
    private boolean        minFlag;
    private boolean        maxFlag;
    private int            limitsMode;
    private double         min;
    private double         max;
    private boolean        lineStem;
    private int            geometryQuality;
    private double         stemWidth;
    private boolean        origOnly;
    private int            glyphType;
}

