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
import llnl.visit.ColorControlPointList;
import llnl.visit.GaussianControlPointList;

// ****************************************************************************
// Class: VolumeAttributes
//
// Purpose:
//    This class contains the plot attributes for the volume plot.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Fri Jan 11 16:25:13 PST 2008
//
// Modifications:
//   
// ****************************************************************************

public class VolumeAttributes extends AttributeSubject implements Plugin
{
    // Enum values
    public final static int RENDERER_SPLATTING = 0;
    public final static int RENDERER_TEXTURE3D = 1;
    public final static int RENDERER_RAYCASTING = 2;
    public final static int RENDERER_RAYCASTINGINTEGRATION = 3;
    public final static int RENDERER_SLIVR = 4;

    public final static int GRADIENTTYPE_CENTEREDDIFFERENCES = 0;
    public final static int GRADIENTTYPE_SOBELOPERATOR = 1;

    public final static int SCALING_LINEAR = 0;
    public final static int SCALING_LOG10 = 1;
    public final static int SCALING_SKEW = 2;

    public final static int SAMPLINGTYPE_KERNELBASED = 0;
    public final static int SAMPLINGTYPE_RASTERIZATION = 1;


    public VolumeAttributes()
    {
        super(26);

        legendFlag = true;
        lightingFlag = true;
        colorControlPoints = new ColorControlPointList();
        opacityAttenuation = 1f;
        freeformFlag = true;
        opacityControlPoints = new GaussianControlPointList();
        resampleTarget = 50000;
        opacityVariable = new String("default");
        freeformOpacity = new byte[256];
        for (int i = 0; i < freeformOpacity.length; ++i)
            freeformOpacity[i] = 0;
        useColorVarMin = false;
        colorVarMin = 0f;
        useColorVarMax = false;
        colorVarMax = 0f;
        useOpacityVarMin = false;
        opacityVarMin = 0f;
        useOpacityVarMax = false;
        opacityVarMax = 0f;
        smoothData = false;
        samplesPerRay = 500;
        rendererType = RENDERER_SPLATTING;
        gradientType = GRADIENTTYPE_SOBELOPERATOR;
        num3DSlices = 200;
        scaling = SCALING_LINEAR;
        skewFactor = 1;
        sampling = SAMPLINGTYPE_RASTERIZATION;
        rendererSamples = 3f;
    }

    public VolumeAttributes(VolumeAttributes obj)
    {
        super(26);

        int i;

        legendFlag = obj.legendFlag;
        lightingFlag = obj.lightingFlag;
        colorControlPoints = new ColorControlPointList(obj.colorControlPoints);
        opacityAttenuation = obj.opacityAttenuation;
        freeformFlag = obj.freeformFlag;
        opacityControlPoints = new GaussianControlPointList(obj.opacityControlPoints);
        resampleTarget = obj.resampleTarget;
        opacityVariable = new String(obj.opacityVariable);
        freeformOpacity = new byte[256];
        for(i = 0; i < obj.freeformOpacity.length; ++i)
            freeformOpacity[i] = obj.freeformOpacity[i];

        useColorVarMin = obj.useColorVarMin;
        colorVarMin = obj.colorVarMin;
        useColorVarMax = obj.useColorVarMax;
        colorVarMax = obj.colorVarMax;
        useOpacityVarMin = obj.useOpacityVarMin;
        opacityVarMin = obj.opacityVarMin;
        useOpacityVarMax = obj.useOpacityVarMax;
        opacityVarMax = obj.opacityVarMax;
        smoothData = obj.smoothData;
        samplesPerRay = obj.samplesPerRay;
        rendererType = obj.rendererType;
        gradientType = obj.gradientType;
        num3DSlices = obj.num3DSlices;
        scaling = obj.scaling;
        skewFactor = obj.skewFactor;
        sampling = obj.sampling;
        rendererSamples = obj.rendererSamples;

        SelectAll();
    }

    public boolean equals(VolumeAttributes obj)
    {
        int i;

        // Compare the freeformOpacity arrays.
        boolean freeformOpacity_equal = true;
        for(i = 0; i < 256 && freeformOpacity_equal; ++i)
            freeformOpacity_equal = (freeformOpacity[i] == obj.freeformOpacity[i]);

        // Create the return value
        return ((legendFlag == obj.legendFlag) &&
                (lightingFlag == obj.lightingFlag) &&
                (colorControlPoints == obj.colorControlPoints) &&
                (opacityAttenuation == obj.opacityAttenuation) &&
                (freeformFlag == obj.freeformFlag) &&
                (opacityControlPoints == obj.opacityControlPoints) &&
                (resampleTarget == obj.resampleTarget) &&
                (opacityVariable == obj.opacityVariable) &&
                freeformOpacity_equal &&
                (useColorVarMin == obj.useColorVarMin) &&
                (colorVarMin == obj.colorVarMin) &&
                (useColorVarMax == obj.useColorVarMax) &&
                (colorVarMax == obj.colorVarMax) &&
                (useOpacityVarMin == obj.useOpacityVarMin) &&
                (opacityVarMin == obj.opacityVarMin) &&
                (useOpacityVarMax == obj.useOpacityVarMax) &&
                (opacityVarMax == obj.opacityVarMax) &&
                (smoothData == obj.smoothData) &&
                (samplesPerRay == obj.samplesPerRay) &&
                (rendererType == obj.rendererType) &&
                (gradientType == obj.gradientType) &&
                (num3DSlices == obj.num3DSlices) &&
                (scaling == obj.scaling) &&
                (skewFactor == obj.skewFactor) &&
                (sampling == obj.sampling) &&
                (rendererSamples == obj.rendererSamples));
    }

    public String GetName() { return "Volume"; }
    public String GetVersion() { return "1.1"; }

    // Property setting methods
    public void SetLegendFlag(boolean legendFlag_)
    {
        legendFlag = legendFlag_;
        Select(0);
    }

    public void SetLightingFlag(boolean lightingFlag_)
    {
        lightingFlag = lightingFlag_;
        Select(1);
    }

    public void SetColorControlPoints(ColorControlPointList colorControlPoints_)
    {
        colorControlPoints = colorControlPoints_;
        Select(2);
    }

    public void SetOpacityAttenuation(float opacityAttenuation_)
    {
        opacityAttenuation = opacityAttenuation_;
        Select(3);
    }

    public void SetFreeformFlag(boolean freeformFlag_)
    {
        freeformFlag = freeformFlag_;
        Select(4);
    }

    public void SetOpacityControlPoints(GaussianControlPointList opacityControlPoints_)
    {
        opacityControlPoints = opacityControlPoints_;
        Select(5);
    }

    public void SetResampleTarget(int resampleTarget_)
    {
        resampleTarget = resampleTarget_;
        Select(6);
    }

    public void SetOpacityVariable(String opacityVariable_)
    {
        opacityVariable = opacityVariable_;
        Select(7);
    }

    public void SetFreeformOpacity(byte[] freeformOpacity_)
    {
        for(int i = 0; i < 256; ++i)
             freeformOpacity[i] = freeformOpacity_[i];
        Select(8);
    }

    public void SetUseColorVarMin(boolean useColorVarMin_)
    {
        useColorVarMin = useColorVarMin_;
        Select(9);
    }

    public void SetColorVarMin(float colorVarMin_)
    {
        colorVarMin = colorVarMin_;
        Select(10);
    }

    public void SetUseColorVarMax(boolean useColorVarMax_)
    {
        useColorVarMax = useColorVarMax_;
        Select(11);
    }

    public void SetColorVarMax(float colorVarMax_)
    {
        colorVarMax = colorVarMax_;
        Select(12);
    }

    public void SetUseOpacityVarMin(boolean useOpacityVarMin_)
    {
        useOpacityVarMin = useOpacityVarMin_;
        Select(13);
    }

    public void SetOpacityVarMin(float opacityVarMin_)
    {
        opacityVarMin = opacityVarMin_;
        Select(14);
    }

    public void SetUseOpacityVarMax(boolean useOpacityVarMax_)
    {
        useOpacityVarMax = useOpacityVarMax_;
        Select(15);
    }

    public void SetOpacityVarMax(float opacityVarMax_)
    {
        opacityVarMax = opacityVarMax_;
        Select(16);
    }

    public void SetSmoothData(boolean smoothData_)
    {
        smoothData = smoothData_;
        Select(17);
    }

    public void SetSamplesPerRay(int samplesPerRay_)
    {
        samplesPerRay = samplesPerRay_;
        Select(18);
    }

    public void SetRendererType(int rendererType_)
    {
        rendererType = rendererType_;
        Select(19);
    }

    public void SetGradientType(int gradientType_)
    {
        gradientType = gradientType_;
        Select(20);
    }

    public void SetNum3DSlices(int num3DSlices_)
    {
        num3DSlices = num3DSlices_;
        Select(21);
    }

    public void SetScaling(int scaling_)
    {
        scaling = scaling_;
        Select(22);
    }

    public void SetSkewFactor(double skewFactor_)
    {
        skewFactor = skewFactor_;
        Select(23);
    }

    public void SetSampling(int sampling_)
    {
        sampling = sampling_;
        Select(24);
    }

    public void SetRendererSamples(float rendererSamples_)
    {
        rendererSamples = rendererSamples_;
        Select(25);
    }

    // Property getting methods
    public boolean                  GetLegendFlag() { return legendFlag; }
    public boolean                  GetLightingFlag() { return lightingFlag; }
    public ColorControlPointList    GetColorControlPoints() { return colorControlPoints; }
    public float                    GetOpacityAttenuation() { return opacityAttenuation; }
    public boolean                  GetFreeformFlag() { return freeformFlag; }
    public GaussianControlPointList GetOpacityControlPoints() { return opacityControlPoints; }
    public int                      GetResampleTarget() { return resampleTarget; }
    public String                   GetOpacityVariable() { return opacityVariable; }
    public byte[]                   GetFreeformOpacity() { return freeformOpacity; }
    public boolean                  GetUseColorVarMin() { return useColorVarMin; }
    public float                    GetColorVarMin() { return colorVarMin; }
    public boolean                  GetUseColorVarMax() { return useColorVarMax; }
    public float                    GetColorVarMax() { return colorVarMax; }
    public boolean                  GetUseOpacityVarMin() { return useOpacityVarMin; }
    public float                    GetOpacityVarMin() { return opacityVarMin; }
    public boolean                  GetUseOpacityVarMax() { return useOpacityVarMax; }
    public float                    GetOpacityVarMax() { return opacityVarMax; }
    public boolean                  GetSmoothData() { return smoothData; }
    public int                      GetSamplesPerRay() { return samplesPerRay; }
    public int                      GetRendererType() { return rendererType; }
    public int                      GetGradientType() { return gradientType; }
    public int                      GetNum3DSlices() { return num3DSlices; }
    public int                      GetScaling() { return scaling; }
    public double                   GetSkewFactor() { return skewFactor; }
    public int                      GetSampling() { return sampling; }
    public float                    GetRendererSamples() { return rendererSamples; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(legendFlag);
        if(WriteSelect(1, buf))
            buf.WriteBool(lightingFlag);
        if(WriteSelect(2, buf))
            colorControlPoints.Write(buf);
        if(WriteSelect(3, buf))
            buf.WriteFloat(opacityAttenuation);
        if(WriteSelect(4, buf))
            buf.WriteBool(freeformFlag);
        if(WriteSelect(5, buf))
            opacityControlPoints.Write(buf);
        if(WriteSelect(6, buf))
            buf.WriteInt(resampleTarget);
        if(WriteSelect(7, buf))
            buf.WriteString(opacityVariable);
        if(WriteSelect(8, buf))
            buf.WriteByteArray(freeformOpacity, true);
        if(WriteSelect(9, buf))
            buf.WriteBool(useColorVarMin);
        if(WriteSelect(10, buf))
            buf.WriteFloat(colorVarMin);
        if(WriteSelect(11, buf))
            buf.WriteBool(useColorVarMax);
        if(WriteSelect(12, buf))
            buf.WriteFloat(colorVarMax);
        if(WriteSelect(13, buf))
            buf.WriteBool(useOpacityVarMin);
        if(WriteSelect(14, buf))
            buf.WriteFloat(opacityVarMin);
        if(WriteSelect(15, buf))
            buf.WriteBool(useOpacityVarMax);
        if(WriteSelect(16, buf))
            buf.WriteFloat(opacityVarMax);
        if(WriteSelect(17, buf))
            buf.WriteBool(smoothData);
        if(WriteSelect(18, buf))
            buf.WriteInt(samplesPerRay);
        if(WriteSelect(19, buf))
            buf.WriteInt(rendererType);
        if(WriteSelect(20, buf))
            buf.WriteInt(gradientType);
        if(WriteSelect(21, buf))
            buf.WriteInt(num3DSlices);
        if(WriteSelect(22, buf))
            buf.WriteInt(scaling);
        if(WriteSelect(23, buf))
            buf.WriteDouble(skewFactor);
        if(WriteSelect(24, buf))
            buf.WriteInt(sampling);
        if(WriteSelect(25, buf))
            buf.WriteFloat(rendererSamples);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetLegendFlag(buf.ReadBool());
                break;
            case 1:
                SetLightingFlag(buf.ReadBool());
                break;
            case 2:
                colorControlPoints.Read(buf);
                Select(2);
                break;
            case 3:
                SetOpacityAttenuation(buf.ReadFloat());
                break;
            case 4:
                SetFreeformFlag(buf.ReadBool());
                break;
            case 5:
                opacityControlPoints.Read(buf);
                Select(5);
                break;
            case 6:
                SetResampleTarget(buf.ReadInt());
                break;
            case 7:
                SetOpacityVariable(buf.ReadString());
                break;
            case 8:
                SetFreeformOpacity(buf.ReadByteArray());
                break;
            case 9:
                SetUseColorVarMin(buf.ReadBool());
                break;
            case 10:
                SetColorVarMin(buf.ReadFloat());
                break;
            case 11:
                SetUseColorVarMax(buf.ReadBool());
                break;
            case 12:
                SetColorVarMax(buf.ReadFloat());
                break;
            case 13:
                SetUseOpacityVarMin(buf.ReadBool());
                break;
            case 14:
                SetOpacityVarMin(buf.ReadFloat());
                break;
            case 15:
                SetUseOpacityVarMax(buf.ReadBool());
                break;
            case 16:
                SetOpacityVarMax(buf.ReadFloat());
                break;
            case 17:
                SetSmoothData(buf.ReadBool());
                break;
            case 18:
                SetSamplesPerRay(buf.ReadInt());
                break;
            case 19:
                SetRendererType(buf.ReadInt());
                break;
            case 20:
                SetGradientType(buf.ReadInt());
                break;
            case 21:
                SetNum3DSlices(buf.ReadInt());
                break;
            case 22:
                SetScaling(buf.ReadInt());
                break;
            case 23:
                SetSkewFactor(buf.ReadDouble());
                break;
            case 24:
                SetSampling(buf.ReadInt());
                break;
            case 25:
                SetRendererSamples(buf.ReadFloat());
                break;
            }
        }
    }


    // Attributes
    private boolean                  legendFlag;
    private boolean                  lightingFlag;
    private ColorControlPointList    colorControlPoints;
    private float                    opacityAttenuation;
    private boolean                  freeformFlag;
    private GaussianControlPointList opacityControlPoints;
    private int                      resampleTarget;
    private String                   opacityVariable;
    private byte[]                   freeformOpacity;
    private boolean                  useColorVarMin;
    private float                    colorVarMin;
    private boolean                  useColorVarMax;
    private float                    colorVarMax;
    private boolean                  useOpacityVarMin;
    private float                    opacityVarMin;
    private boolean                  useOpacityVarMax;
    private float                    opacityVarMax;
    private boolean                  smoothData;
    private int                      samplesPerRay;
    private int                      rendererType;
    private int                      gradientType;
    private int                      num3DSlices;
    private int                      scaling;
    private double                   skewFactor;
    private int                      sampling;
    private float                    rendererSamples;
}

