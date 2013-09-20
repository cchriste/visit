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

package llnl.visit.plots;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: PseudocolorAttributes
//
// Purpose:
//    Attributes for the pseudocolor plot
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class PseudocolorAttributes extends AttributeSubject implements Plugin
{
    private static int PseudocolorAttributes_numAdditionalAtts = 44;

    // Enum values
    public final static int SCALING_LINEAR = 0;
    public final static int SCALING_LOG = 1;
    public final static int SCALING_SKEW = 2;

    public final static int LIMITSMODE_ORIGINALDATA = 0;
    public final static int LIMITSMODE_CURRENTPLOT = 1;

    public final static int CENTERING_NATURAL = 0;
    public final static int CENTERING_NODAL = 1;
    public final static int CENTERING_ZONAL = 2;

    public final static int OPACITYTYPE_COLORTABLE = 0;
    public final static int OPACITYTYPE_FULLYOPAQUE = 1;
    public final static int OPACITYTYPE_CONSTANT = 2;
    public final static int OPACITYTYPE_RAMP = 3;
    public final static int OPACITYTYPE_VARIABLERANGE = 4;

    public final static int POINTTYPE_BOX = 0;
    public final static int POINTTYPE_AXIS = 1;
    public final static int POINTTYPE_ICOSAHEDRON = 2;
    public final static int POINTTYPE_OCTAHEDRON = 3;
    public final static int POINTTYPE_TETRAHEDRON = 4;
    public final static int POINTTYPE_SPHEREGEOMETRY = 5;
    public final static int POINTTYPE_POINT = 6;
    public final static int POINTTYPE_SPHERE = 7;

    public final static int LINETYPE_LINE = 0;
    public final static int LINETYPE_TUBE = 1;
    public final static int LINETYPE_RIBBON = 2;

    public final static int ENDPOINTTYPE_NONE = 0;
    public final static int ENDPOINTTYPE_TAILS = 1;
    public final static int ENDPOINTTYPE_HEADS = 2;
    public final static int ENDPOINTTYPE_BOTH = 3;

    public final static int ENDPOINTSTYLE_SPHERES = 0;
    public final static int ENDPOINTSTYLE_CONES = 1;

    public final static int SIZETYPE_ABSOLUTE = 0;
    public final static int SIZETYPE_FRACTIONOFBBOX = 1;


    public PseudocolorAttributes()
    {
        super(PseudocolorAttributes_numAdditionalAtts);

        scaling = SCALING_LINEAR;
        skewFactor = 1;
        limitsMode = LIMITSMODE_ORIGINALDATA;
        minFlag = false;
        min = 0;
        maxFlag = false;
        max = 1;
        centering = CENTERING_NATURAL;
        colorTableName = new String("hot");
        invertColorTable = false;
        opacityType = OPACITYTYPE_FULLYOPAQUE;
        opacityVariable = new String("");
        opacity = 1;
        opacityVarMin = 0;
        opacityVarMax = 1;
        opacityVarMinFlag = false;
        opacityVarMaxFlag = false;
        pointSize = 0.05;
        pointType = POINTTYPE_POINT;
        pointSizeVarEnabled = false;
        pointSizeVar = new String("default");
        pointSizePixels = 2;
        lineType = LINETYPE_LINE;
        lineStyle = 0;
        lineWidth = 0;
        tubeDisplayDensity = 10;
        tubeRadiusSizeType = SIZETYPE_FRACTIONOFBBOX;
        tubeRadiusAbsolute = 0.125;
        tubeRadiusBBox = 0.005;
        varyTubeRadius = false;
        varyTubeRadiusVariable = new String("");
        varyTubeRadiusFactor = 10;
        endPointType = ENDPOINTTYPE_NONE;
        endPointStyle = ENDPOINTSTYLE_SPHERES;
        endPointRadiusSizeType = SIZETYPE_FRACTIONOFBBOX;
        endPointRadiusAbsolute = 1;
        endPointRadiusBBox = 0.005;
        endPointRatio = 2;
        renderSurfaces = 1;
        renderWireframe = 0;
        renderPoints = 0;
        smoothingLevel = 0;
        legendFlag = true;
        lightingFlag = true;
    }

    public PseudocolorAttributes(int nMoreFields)
    {
        super(PseudocolorAttributes_numAdditionalAtts + nMoreFields);

        scaling = SCALING_LINEAR;
        skewFactor = 1;
        limitsMode = LIMITSMODE_ORIGINALDATA;
        minFlag = false;
        min = 0;
        maxFlag = false;
        max = 1;
        centering = CENTERING_NATURAL;
        colorTableName = new String("hot");
        invertColorTable = false;
        opacityType = OPACITYTYPE_FULLYOPAQUE;
        opacityVariable = new String("");
        opacity = 1;
        opacityVarMin = 0;
        opacityVarMax = 1;
        opacityVarMinFlag = false;
        opacityVarMaxFlag = false;
        pointSize = 0.05;
        pointType = POINTTYPE_POINT;
        pointSizeVarEnabled = false;
        pointSizeVar = new String("default");
        pointSizePixels = 2;
        lineType = LINETYPE_LINE;
        lineStyle = 0;
        lineWidth = 0;
        tubeDisplayDensity = 10;
        tubeRadiusSizeType = SIZETYPE_FRACTIONOFBBOX;
        tubeRadiusAbsolute = 0.125;
        tubeRadiusBBox = 0.005;
        varyTubeRadius = false;
        varyTubeRadiusVariable = new String("");
        varyTubeRadiusFactor = 10;
        endPointType = ENDPOINTTYPE_NONE;
        endPointStyle = ENDPOINTSTYLE_SPHERES;
        endPointRadiusSizeType = SIZETYPE_FRACTIONOFBBOX;
        endPointRadiusAbsolute = 1;
        endPointRadiusBBox = 0.005;
        endPointRatio = 2;
        renderSurfaces = 1;
        renderWireframe = 0;
        renderPoints = 0;
        smoothingLevel = 0;
        legendFlag = true;
        lightingFlag = true;
    }

    public PseudocolorAttributes(PseudocolorAttributes obj)
    {
        super(PseudocolorAttributes_numAdditionalAtts);

        scaling = obj.scaling;
        skewFactor = obj.skewFactor;
        limitsMode = obj.limitsMode;
        minFlag = obj.minFlag;
        min = obj.min;
        maxFlag = obj.maxFlag;
        max = obj.max;
        centering = obj.centering;
        colorTableName = new String(obj.colorTableName);
        invertColorTable = obj.invertColorTable;
        opacityType = obj.opacityType;
        opacityVariable = new String(obj.opacityVariable);
        opacity = obj.opacity;
        opacityVarMin = obj.opacityVarMin;
        opacityVarMax = obj.opacityVarMax;
        opacityVarMinFlag = obj.opacityVarMinFlag;
        opacityVarMaxFlag = obj.opacityVarMaxFlag;
        pointSize = obj.pointSize;
        pointType = obj.pointType;
        pointSizeVarEnabled = obj.pointSizeVarEnabled;
        pointSizeVar = new String(obj.pointSizeVar);
        pointSizePixels = obj.pointSizePixels;
        lineType = obj.lineType;
        lineStyle = obj.lineStyle;
        lineWidth = obj.lineWidth;
        tubeDisplayDensity = obj.tubeDisplayDensity;
        tubeRadiusSizeType = obj.tubeRadiusSizeType;
        tubeRadiusAbsolute = obj.tubeRadiusAbsolute;
        tubeRadiusBBox = obj.tubeRadiusBBox;
        varyTubeRadius = obj.varyTubeRadius;
        varyTubeRadiusVariable = new String(obj.varyTubeRadiusVariable);
        varyTubeRadiusFactor = obj.varyTubeRadiusFactor;
        endPointType = obj.endPointType;
        endPointStyle = obj.endPointStyle;
        endPointRadiusSizeType = obj.endPointRadiusSizeType;
        endPointRadiusAbsolute = obj.endPointRadiusAbsolute;
        endPointRadiusBBox = obj.endPointRadiusBBox;
        endPointRatio = obj.endPointRatio;
        renderSurfaces = obj.renderSurfaces;
        renderWireframe = obj.renderWireframe;
        renderPoints = obj.renderPoints;
        smoothingLevel = obj.smoothingLevel;
        legendFlag = obj.legendFlag;
        lightingFlag = obj.lightingFlag;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return PseudocolorAttributes_numAdditionalAtts;
    }

    public boolean equals(PseudocolorAttributes obj)
    {
        // Create the return value
        return ((scaling == obj.scaling) &&
                (skewFactor == obj.skewFactor) &&
                (limitsMode == obj.limitsMode) &&
                (minFlag == obj.minFlag) &&
                (min == obj.min) &&
                (maxFlag == obj.maxFlag) &&
                (max == obj.max) &&
                (centering == obj.centering) &&
                (colorTableName.equals(obj.colorTableName)) &&
                (invertColorTable == obj.invertColorTable) &&
                (opacityType == obj.opacityType) &&
                (opacityVariable.equals(obj.opacityVariable)) &&
                (opacity == obj.opacity) &&
                (opacityVarMin == obj.opacityVarMin) &&
                (opacityVarMax == obj.opacityVarMax) &&
                (opacityVarMinFlag == obj.opacityVarMinFlag) &&
                (opacityVarMaxFlag == obj.opacityVarMaxFlag) &&
                (pointSize == obj.pointSize) &&
                (pointType == obj.pointType) &&
                (pointSizeVarEnabled == obj.pointSizeVarEnabled) &&
                (pointSizeVar.equals(obj.pointSizeVar)) &&
                (pointSizePixels == obj.pointSizePixels) &&
                (lineType == obj.lineType) &&
                (lineStyle == obj.lineStyle) &&
                (lineWidth == obj.lineWidth) &&
                (tubeDisplayDensity == obj.tubeDisplayDensity) &&
                (tubeRadiusSizeType == obj.tubeRadiusSizeType) &&
                (tubeRadiusAbsolute == obj.tubeRadiusAbsolute) &&
                (tubeRadiusBBox == obj.tubeRadiusBBox) &&
                (varyTubeRadius == obj.varyTubeRadius) &&
                (varyTubeRadiusVariable.equals(obj.varyTubeRadiusVariable)) &&
                (varyTubeRadiusFactor == obj.varyTubeRadiusFactor) &&
                (endPointType == obj.endPointType) &&
                (endPointStyle == obj.endPointStyle) &&
                (endPointRadiusSizeType == obj.endPointRadiusSizeType) &&
                (endPointRadiusAbsolute == obj.endPointRadiusAbsolute) &&
                (endPointRadiusBBox == obj.endPointRadiusBBox) &&
                (endPointRatio == obj.endPointRatio) &&
                (renderSurfaces == obj.renderSurfaces) &&
                (renderWireframe == obj.renderWireframe) &&
                (renderPoints == obj.renderPoints) &&
                (smoothingLevel == obj.smoothingLevel) &&
                (legendFlag == obj.legendFlag) &&
                (lightingFlag == obj.lightingFlag));
    }

    public String GetName() { return "Pseudocolor"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetScaling(int scaling_)
    {
        scaling = scaling_;
        Select(0);
    }

    public void SetSkewFactor(double skewFactor_)
    {
        skewFactor = skewFactor_;
        Select(1);
    }

    public void SetLimitsMode(int limitsMode_)
    {
        limitsMode = limitsMode_;
        Select(2);
    }

    public void SetMinFlag(boolean minFlag_)
    {
        minFlag = minFlag_;
        Select(3);
    }

    public void SetMin(double min_)
    {
        min = min_;
        Select(4);
    }

    public void SetMaxFlag(boolean maxFlag_)
    {
        maxFlag = maxFlag_;
        Select(5);
    }

    public void SetMax(double max_)
    {
        max = max_;
        Select(6);
    }

    public void SetCentering(int centering_)
    {
        centering = centering_;
        Select(7);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(8);
    }

    public void SetInvertColorTable(boolean invertColorTable_)
    {
        invertColorTable = invertColorTable_;
        Select(9);
    }

    public void SetOpacityType(int opacityType_)
    {
        opacityType = opacityType_;
        Select(10);
    }

    public void SetOpacityVariable(String opacityVariable_)
    {
        opacityVariable = opacityVariable_;
        Select(11);
    }

    public void SetOpacity(double opacity_)
    {
        opacity = opacity_;
        Select(12);
    }

    public void SetOpacityVarMin(double opacityVarMin_)
    {
        opacityVarMin = opacityVarMin_;
        Select(13);
    }

    public void SetOpacityVarMax(double opacityVarMax_)
    {
        opacityVarMax = opacityVarMax_;
        Select(14);
    }

    public void SetOpacityVarMinFlag(boolean opacityVarMinFlag_)
    {
        opacityVarMinFlag = opacityVarMinFlag_;
        Select(15);
    }

    public void SetOpacityVarMaxFlag(boolean opacityVarMaxFlag_)
    {
        opacityVarMaxFlag = opacityVarMaxFlag_;
        Select(16);
    }

    public void SetPointSize(double pointSize_)
    {
        pointSize = pointSize_;
        Select(17);
    }

    public void SetPointType(int pointType_)
    {
        pointType = pointType_;
        Select(18);
    }

    public void SetPointSizeVarEnabled(boolean pointSizeVarEnabled_)
    {
        pointSizeVarEnabled = pointSizeVarEnabled_;
        Select(19);
    }

    public void SetPointSizeVar(String pointSizeVar_)
    {
        pointSizeVar = pointSizeVar_;
        Select(20);
    }

    public void SetPointSizePixels(int pointSizePixels_)
    {
        pointSizePixels = pointSizePixels_;
        Select(21);
    }

    public void SetLineType(int lineType_)
    {
        lineType = lineType_;
        Select(22);
    }

    public void SetLineStyle(int lineStyle_)
    {
        lineStyle = lineStyle_;
        Select(23);
    }

    public void SetLineWidth(int lineWidth_)
    {
        lineWidth = lineWidth_;
        Select(24);
    }

    public void SetTubeDisplayDensity(int tubeDisplayDensity_)
    {
        tubeDisplayDensity = tubeDisplayDensity_;
        Select(25);
    }

    public void SetTubeRadiusSizeType(int tubeRadiusSizeType_)
    {
        tubeRadiusSizeType = tubeRadiusSizeType_;
        Select(26);
    }

    public void SetTubeRadiusAbsolute(double tubeRadiusAbsolute_)
    {
        tubeRadiusAbsolute = tubeRadiusAbsolute_;
        Select(27);
    }

    public void SetTubeRadiusBBox(double tubeRadiusBBox_)
    {
        tubeRadiusBBox = tubeRadiusBBox_;
        Select(28);
    }

    public void SetVaryTubeRadius(boolean varyTubeRadius_)
    {
        varyTubeRadius = varyTubeRadius_;
        Select(29);
    }

    public void SetVaryTubeRadiusVariable(String varyTubeRadiusVariable_)
    {
        varyTubeRadiusVariable = varyTubeRadiusVariable_;
        Select(30);
    }

    public void SetVaryTubeRadiusFactor(double varyTubeRadiusFactor_)
    {
        varyTubeRadiusFactor = varyTubeRadiusFactor_;
        Select(31);
    }

    public void SetEndPointType(int endPointType_)
    {
        endPointType = endPointType_;
        Select(32);
    }

    public void SetEndPointStyle(int endPointStyle_)
    {
        endPointStyle = endPointStyle_;
        Select(33);
    }

    public void SetEndPointRadiusSizeType(int endPointRadiusSizeType_)
    {
        endPointRadiusSizeType = endPointRadiusSizeType_;
        Select(34);
    }

    public void SetEndPointRadiusAbsolute(double endPointRadiusAbsolute_)
    {
        endPointRadiusAbsolute = endPointRadiusAbsolute_;
        Select(35);
    }

    public void SetEndPointRadiusBBox(double endPointRadiusBBox_)
    {
        endPointRadiusBBox = endPointRadiusBBox_;
        Select(36);
    }

    public void SetEndPointRatio(double endPointRatio_)
    {
        endPointRatio = endPointRatio_;
        Select(37);
    }

    public void SetRenderSurfaces(int renderSurfaces_)
    {
        renderSurfaces = renderSurfaces_;
        Select(38);
    }

    public void SetRenderWireframe(int renderWireframe_)
    {
        renderWireframe = renderWireframe_;
        Select(39);
    }

    public void SetRenderPoints(int renderPoints_)
    {
        renderPoints = renderPoints_;
        Select(40);
    }

    public void SetSmoothingLevel(int smoothingLevel_)
    {
        smoothingLevel = smoothingLevel_;
        Select(41);
    }

    public void SetLegendFlag(boolean legendFlag_)
    {
        legendFlag = legendFlag_;
        Select(42);
    }

    public void SetLightingFlag(boolean lightingFlag_)
    {
        lightingFlag = lightingFlag_;
        Select(43);
    }

    // Property getting methods
    public int     GetScaling() { return scaling; }
    public double  GetSkewFactor() { return skewFactor; }
    public int     GetLimitsMode() { return limitsMode; }
    public boolean GetMinFlag() { return minFlag; }
    public double  GetMin() { return min; }
    public boolean GetMaxFlag() { return maxFlag; }
    public double  GetMax() { return max; }
    public int     GetCentering() { return centering; }
    public String  GetColorTableName() { return colorTableName; }
    public boolean GetInvertColorTable() { return invertColorTable; }
    public int     GetOpacityType() { return opacityType; }
    public String  GetOpacityVariable() { return opacityVariable; }
    public double  GetOpacity() { return opacity; }
    public double  GetOpacityVarMin() { return opacityVarMin; }
    public double  GetOpacityVarMax() { return opacityVarMax; }
    public boolean GetOpacityVarMinFlag() { return opacityVarMinFlag; }
    public boolean GetOpacityVarMaxFlag() { return opacityVarMaxFlag; }
    public double  GetPointSize() { return pointSize; }
    public int     GetPointType() { return pointType; }
    public boolean GetPointSizeVarEnabled() { return pointSizeVarEnabled; }
    public String  GetPointSizeVar() { return pointSizeVar; }
    public int     GetPointSizePixels() { return pointSizePixels; }
    public int     GetLineType() { return lineType; }
    public int     GetLineStyle() { return lineStyle; }
    public int     GetLineWidth() { return lineWidth; }
    public int     GetTubeDisplayDensity() { return tubeDisplayDensity; }
    public int     GetTubeRadiusSizeType() { return tubeRadiusSizeType; }
    public double  GetTubeRadiusAbsolute() { return tubeRadiusAbsolute; }
    public double  GetTubeRadiusBBox() { return tubeRadiusBBox; }
    public boolean GetVaryTubeRadius() { return varyTubeRadius; }
    public String  GetVaryTubeRadiusVariable() { return varyTubeRadiusVariable; }
    public double  GetVaryTubeRadiusFactor() { return varyTubeRadiusFactor; }
    public int     GetEndPointType() { return endPointType; }
    public int     GetEndPointStyle() { return endPointStyle; }
    public int     GetEndPointRadiusSizeType() { return endPointRadiusSizeType; }
    public double  GetEndPointRadiusAbsolute() { return endPointRadiusAbsolute; }
    public double  GetEndPointRadiusBBox() { return endPointRadiusBBox; }
    public double  GetEndPointRatio() { return endPointRatio; }
    public int     GetRenderSurfaces() { return renderSurfaces; }
    public int     GetRenderWireframe() { return renderWireframe; }
    public int     GetRenderPoints() { return renderPoints; }
    public int     GetSmoothingLevel() { return smoothingLevel; }
    public boolean GetLegendFlag() { return legendFlag; }
    public boolean GetLightingFlag() { return lightingFlag; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(scaling);
        if(WriteSelect(1, buf))
            buf.WriteDouble(skewFactor);
        if(WriteSelect(2, buf))
            buf.WriteInt(limitsMode);
        if(WriteSelect(3, buf))
            buf.WriteBool(minFlag);
        if(WriteSelect(4, buf))
            buf.WriteDouble(min);
        if(WriteSelect(5, buf))
            buf.WriteBool(maxFlag);
        if(WriteSelect(6, buf))
            buf.WriteDouble(max);
        if(WriteSelect(7, buf))
            buf.WriteInt(centering);
        if(WriteSelect(8, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(9, buf))
            buf.WriteBool(invertColorTable);
        if(WriteSelect(10, buf))
            buf.WriteInt(opacityType);
        if(WriteSelect(11, buf))
            buf.WriteString(opacityVariable);
        if(WriteSelect(12, buf))
            buf.WriteDouble(opacity);
        if(WriteSelect(13, buf))
            buf.WriteDouble(opacityVarMin);
        if(WriteSelect(14, buf))
            buf.WriteDouble(opacityVarMax);
        if(WriteSelect(15, buf))
            buf.WriteBool(opacityVarMinFlag);
        if(WriteSelect(16, buf))
            buf.WriteBool(opacityVarMaxFlag);
        if(WriteSelect(17, buf))
            buf.WriteDouble(pointSize);
        if(WriteSelect(18, buf))
            buf.WriteInt(pointType);
        if(WriteSelect(19, buf))
            buf.WriteBool(pointSizeVarEnabled);
        if(WriteSelect(20, buf))
            buf.WriteString(pointSizeVar);
        if(WriteSelect(21, buf))
            buf.WriteInt(pointSizePixels);
        if(WriteSelect(22, buf))
            buf.WriteInt(lineType);
        if(WriteSelect(23, buf))
            buf.WriteInt(lineStyle);
        if(WriteSelect(24, buf))
            buf.WriteInt(lineWidth);
        if(WriteSelect(25, buf))
            buf.WriteInt(tubeDisplayDensity);
        if(WriteSelect(26, buf))
            buf.WriteInt(tubeRadiusSizeType);
        if(WriteSelect(27, buf))
            buf.WriteDouble(tubeRadiusAbsolute);
        if(WriteSelect(28, buf))
            buf.WriteDouble(tubeRadiusBBox);
        if(WriteSelect(29, buf))
            buf.WriteBool(varyTubeRadius);
        if(WriteSelect(30, buf))
            buf.WriteString(varyTubeRadiusVariable);
        if(WriteSelect(31, buf))
            buf.WriteDouble(varyTubeRadiusFactor);
        if(WriteSelect(32, buf))
            buf.WriteInt(endPointType);
        if(WriteSelect(33, buf))
            buf.WriteInt(endPointStyle);
        if(WriteSelect(34, buf))
            buf.WriteInt(endPointRadiusSizeType);
        if(WriteSelect(35, buf))
            buf.WriteDouble(endPointRadiusAbsolute);
        if(WriteSelect(36, buf))
            buf.WriteDouble(endPointRadiusBBox);
        if(WriteSelect(37, buf))
            buf.WriteDouble(endPointRatio);
        if(WriteSelect(38, buf))
            buf.WriteInt(renderSurfaces);
        if(WriteSelect(39, buf))
            buf.WriteInt(renderWireframe);
        if(WriteSelect(40, buf))
            buf.WriteInt(renderPoints);
        if(WriteSelect(41, buf))
            buf.WriteInt(smoothingLevel);
        if(WriteSelect(42, buf))
            buf.WriteBool(legendFlag);
        if(WriteSelect(43, buf))
            buf.WriteBool(lightingFlag);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetScaling(buf.ReadInt());
            break;
        case 1:
            SetSkewFactor(buf.ReadDouble());
            break;
        case 2:
            SetLimitsMode(buf.ReadInt());
            break;
        case 3:
            SetMinFlag(buf.ReadBool());
            break;
        case 4:
            SetMin(buf.ReadDouble());
            break;
        case 5:
            SetMaxFlag(buf.ReadBool());
            break;
        case 6:
            SetMax(buf.ReadDouble());
            break;
        case 7:
            SetCentering(buf.ReadInt());
            break;
        case 8:
            SetColorTableName(buf.ReadString());
            break;
        case 9:
            SetInvertColorTable(buf.ReadBool());
            break;
        case 10:
            SetOpacityType(buf.ReadInt());
            break;
        case 11:
            SetOpacityVariable(buf.ReadString());
            break;
        case 12:
            SetOpacity(buf.ReadDouble());
            break;
        case 13:
            SetOpacityVarMin(buf.ReadDouble());
            break;
        case 14:
            SetOpacityVarMax(buf.ReadDouble());
            break;
        case 15:
            SetOpacityVarMinFlag(buf.ReadBool());
            break;
        case 16:
            SetOpacityVarMaxFlag(buf.ReadBool());
            break;
        case 17:
            SetPointSize(buf.ReadDouble());
            break;
        case 18:
            SetPointType(buf.ReadInt());
            break;
        case 19:
            SetPointSizeVarEnabled(buf.ReadBool());
            break;
        case 20:
            SetPointSizeVar(buf.ReadString());
            break;
        case 21:
            SetPointSizePixels(buf.ReadInt());
            break;
        case 22:
            SetLineType(buf.ReadInt());
            break;
        case 23:
            SetLineStyle(buf.ReadInt());
            break;
        case 24:
            SetLineWidth(buf.ReadInt());
            break;
        case 25:
            SetTubeDisplayDensity(buf.ReadInt());
            break;
        case 26:
            SetTubeRadiusSizeType(buf.ReadInt());
            break;
        case 27:
            SetTubeRadiusAbsolute(buf.ReadDouble());
            break;
        case 28:
            SetTubeRadiusBBox(buf.ReadDouble());
            break;
        case 29:
            SetVaryTubeRadius(buf.ReadBool());
            break;
        case 30:
            SetVaryTubeRadiusVariable(buf.ReadString());
            break;
        case 31:
            SetVaryTubeRadiusFactor(buf.ReadDouble());
            break;
        case 32:
            SetEndPointType(buf.ReadInt());
            break;
        case 33:
            SetEndPointStyle(buf.ReadInt());
            break;
        case 34:
            SetEndPointRadiusSizeType(buf.ReadInt());
            break;
        case 35:
            SetEndPointRadiusAbsolute(buf.ReadDouble());
            break;
        case 36:
            SetEndPointRadiusBBox(buf.ReadDouble());
            break;
        case 37:
            SetEndPointRatio(buf.ReadDouble());
            break;
        case 38:
            SetRenderSurfaces(buf.ReadInt());
            break;
        case 39:
            SetRenderWireframe(buf.ReadInt());
            break;
        case 40:
            SetRenderPoints(buf.ReadInt());
            break;
        case 41:
            SetSmoothingLevel(buf.ReadInt());
            break;
        case 42:
            SetLegendFlag(buf.ReadBool());
            break;
        case 43:
            SetLightingFlag(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "scaling = ";
        if(scaling == SCALING_LINEAR)
            str = str + "SCALING_LINEAR";
        if(scaling == SCALING_LOG)
            str = str + "SCALING_LOG";
        if(scaling == SCALING_SKEW)
            str = str + "SCALING_SKEW";
        str = str + "\n";
        str = str + doubleToString("skewFactor", skewFactor, indent) + "\n";
        str = str + indent + "limitsMode = ";
        if(limitsMode == LIMITSMODE_ORIGINALDATA)
            str = str + "LIMITSMODE_ORIGINALDATA";
        if(limitsMode == LIMITSMODE_CURRENTPLOT)
            str = str + "LIMITSMODE_CURRENTPLOT";
        str = str + "\n";
        str = str + boolToString("minFlag", minFlag, indent) + "\n";
        str = str + doubleToString("min", min, indent) + "\n";
        str = str + boolToString("maxFlag", maxFlag, indent) + "\n";
        str = str + doubleToString("max", max, indent) + "\n";
        str = str + indent + "centering = ";
        if(centering == CENTERING_NATURAL)
            str = str + "CENTERING_NATURAL";
        if(centering == CENTERING_NODAL)
            str = str + "CENTERING_NODAL";
        if(centering == CENTERING_ZONAL)
            str = str + "CENTERING_ZONAL";
        str = str + "\n";
        str = str + stringToString("colorTableName", colorTableName, indent) + "\n";
        str = str + boolToString("invertColorTable", invertColorTable, indent) + "\n";
        str = str + indent + "opacityType = ";
        if(opacityType == OPACITYTYPE_COLORTABLE)
            str = str + "OPACITYTYPE_COLORTABLE";
        if(opacityType == OPACITYTYPE_FULLYOPAQUE)
            str = str + "OPACITYTYPE_FULLYOPAQUE";
        if(opacityType == OPACITYTYPE_CONSTANT)
            str = str + "OPACITYTYPE_CONSTANT";
        if(opacityType == OPACITYTYPE_RAMP)
            str = str + "OPACITYTYPE_RAMP";
        if(opacityType == OPACITYTYPE_VARIABLERANGE)
            str = str + "OPACITYTYPE_VARIABLERANGE";
        str = str + "\n";
        str = str + stringToString("opacityVariable", opacityVariable, indent) + "\n";
        str = str + doubleToString("opacity", opacity, indent) + "\n";
        str = str + doubleToString("opacityVarMin", opacityVarMin, indent) + "\n";
        str = str + doubleToString("opacityVarMax", opacityVarMax, indent) + "\n";
        str = str + boolToString("opacityVarMinFlag", opacityVarMinFlag, indent) + "\n";
        str = str + boolToString("opacityVarMaxFlag", opacityVarMaxFlag, indent) + "\n";
        str = str + doubleToString("pointSize", pointSize, indent) + "\n";
        str = str + indent + "pointType = ";
        if(pointType == POINTTYPE_BOX)
            str = str + "POINTTYPE_BOX";
        if(pointType == POINTTYPE_AXIS)
            str = str + "POINTTYPE_AXIS";
        if(pointType == POINTTYPE_ICOSAHEDRON)
            str = str + "POINTTYPE_ICOSAHEDRON";
        if(pointType == POINTTYPE_OCTAHEDRON)
            str = str + "POINTTYPE_OCTAHEDRON";
        if(pointType == POINTTYPE_TETRAHEDRON)
            str = str + "POINTTYPE_TETRAHEDRON";
        if(pointType == POINTTYPE_SPHEREGEOMETRY)
            str = str + "POINTTYPE_SPHEREGEOMETRY";
        if(pointType == POINTTYPE_POINT)
            str = str + "POINTTYPE_POINT";
        if(pointType == POINTTYPE_SPHERE)
            str = str + "POINTTYPE_SPHERE";
        str = str + "\n";
        str = str + boolToString("pointSizeVarEnabled", pointSizeVarEnabled, indent) + "\n";
        str = str + stringToString("pointSizeVar", pointSizeVar, indent) + "\n";
        str = str + intToString("pointSizePixels", pointSizePixels, indent) + "\n";
        str = str + indent + "lineType = ";
        if(lineType == LINETYPE_LINE)
            str = str + "LINETYPE_LINE";
        if(lineType == LINETYPE_TUBE)
            str = str + "LINETYPE_TUBE";
        if(lineType == LINETYPE_RIBBON)
            str = str + "LINETYPE_RIBBON";
        str = str + "\n";
        str = str + intToString("lineStyle", lineStyle, indent) + "\n";
        str = str + intToString("lineWidth", lineWidth, indent) + "\n";
        str = str + intToString("tubeDisplayDensity", tubeDisplayDensity, indent) + "\n";
        str = str + indent + "tubeRadiusSizeType = ";
        if(tubeRadiusSizeType == SIZETYPE_ABSOLUTE)
            str = str + "SIZETYPE_ABSOLUTE";
        if(tubeRadiusSizeType == SIZETYPE_FRACTIONOFBBOX)
            str = str + "SIZETYPE_FRACTIONOFBBOX";
        str = str + "\n";
        str = str + doubleToString("tubeRadiusAbsolute", tubeRadiusAbsolute, indent) + "\n";
        str = str + doubleToString("tubeRadiusBBox", tubeRadiusBBox, indent) + "\n";
        str = str + boolToString("varyTubeRadius", varyTubeRadius, indent) + "\n";
        str = str + stringToString("varyTubeRadiusVariable", varyTubeRadiusVariable, indent) + "\n";
        str = str + doubleToString("varyTubeRadiusFactor", varyTubeRadiusFactor, indent) + "\n";
        str = str + indent + "endPointType = ";
        if(endPointType == ENDPOINTTYPE_NONE)
            str = str + "ENDPOINTTYPE_NONE";
        if(endPointType == ENDPOINTTYPE_TAILS)
            str = str + "ENDPOINTTYPE_TAILS";
        if(endPointType == ENDPOINTTYPE_HEADS)
            str = str + "ENDPOINTTYPE_HEADS";
        if(endPointType == ENDPOINTTYPE_BOTH)
            str = str + "ENDPOINTTYPE_BOTH";
        str = str + "\n";
        str = str + indent + "endPointStyle = ";
        if(endPointStyle == ENDPOINTSTYLE_SPHERES)
            str = str + "ENDPOINTSTYLE_SPHERES";
        if(endPointStyle == ENDPOINTSTYLE_CONES)
            str = str + "ENDPOINTSTYLE_CONES";
        str = str + "\n";
        str = str + indent + "endPointRadiusSizeType = ";
        if(endPointRadiusSizeType == SIZETYPE_ABSOLUTE)
            str = str + "SIZETYPE_ABSOLUTE";
        if(endPointRadiusSizeType == SIZETYPE_FRACTIONOFBBOX)
            str = str + "SIZETYPE_FRACTIONOFBBOX";
        str = str + "\n";
        str = str + doubleToString("endPointRadiusAbsolute", endPointRadiusAbsolute, indent) + "\n";
        str = str + doubleToString("endPointRadiusBBox", endPointRadiusBBox, indent) + "\n";
        str = str + doubleToString("endPointRatio", endPointRatio, indent) + "\n";
        str = str + intToString("renderSurfaces", renderSurfaces, indent) + "\n";
        str = str + intToString("renderWireframe", renderWireframe, indent) + "\n";
        str = str + intToString("renderPoints", renderPoints, indent) + "\n";
        str = str + intToString("smoothingLevel", smoothingLevel, indent) + "\n";
        str = str + boolToString("legendFlag", legendFlag, indent) + "\n";
        str = str + boolToString("lightingFlag", lightingFlag, indent) + "\n";
        return str;
    }


    // Attributes
    private int     scaling;
    private double  skewFactor;
    private int     limitsMode;
    private boolean minFlag;
    private double  min;
    private boolean maxFlag;
    private double  max;
    private int     centering;
    private String  colorTableName;
    private boolean invertColorTable;
    private int     opacityType;
    private String  opacityVariable;
    private double  opacity;
    private double  opacityVarMin;
    private double  opacityVarMax;
    private boolean opacityVarMinFlag;
    private boolean opacityVarMaxFlag;
    private double  pointSize;
    private int     pointType;
    private boolean pointSizeVarEnabled;
    private String  pointSizeVar;
    private int     pointSizePixels;
    private int     lineType;
    private int     lineStyle;
    private int     lineWidth;
    private int     tubeDisplayDensity;
    private int     tubeRadiusSizeType;
    private double  tubeRadiusAbsolute;
    private double  tubeRadiusBBox;
    private boolean varyTubeRadius;
    private String  varyTubeRadiusVariable;
    private double  varyTubeRadiusFactor;
    private int     endPointType;
    private int     endPointStyle;
    private int     endPointRadiusSizeType;
    private double  endPointRadiusAbsolute;
    private double  endPointRadiusBBox;
    private double  endPointRatio;
    private int     renderSurfaces;
    private int     renderWireframe;
    private int     renderPoints;
    private int     smoothingLevel;
    private boolean legendFlag;
    private boolean lightingFlag;
}

