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
// Class: VectorAttributes
//
// Purpose:
//    Attributes for the vector plot
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Mon Mar 19 16:19:48 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class VectorAttributes extends AttributeSubject implements Plugin
{
    // Enum values
    public final static int ORIGINTYPE_HEAD = 0;
    public final static int ORIGINTYPE_MIDDLE = 1;
    public final static int ORIGINTYPE_TAIL = 2;

    public final static int LIMITSMODE_ORIGINALDATA = 0;
    public final static int LIMITSMODE_CURRENTPLOT = 1;


    public VectorAttributes()
    {
        super(23);

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
        vectorOrigin = ORIGINTYPE_TAIL;
        minFlag = false;
        maxFlag = false;
        limitsMode = LIMITSMODE_ORIGINALDATA;
        min = 0;
        max = 1;
        lineStem = true;
        highQuality = false;
        stemWidth = 0.08;
    }

    public VectorAttributes(VectorAttributes obj)
    {
        super(23);

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
        vectorOrigin = obj.vectorOrigin;
        minFlag = obj.minFlag;
        maxFlag = obj.maxFlag;
        limitsMode = obj.limitsMode;
        min = obj.min;
        max = obj.max;
        lineStem = obj.lineStem;
        highQuality = obj.highQuality;
        stemWidth = obj.stemWidth;

        SelectAll();
    }

    public boolean equals(VectorAttributes obj)
    {
        // Create the return value
        return ((useStride == obj.useStride) &&
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
                (colorTableName == obj.colorTableName) &&
                (vectorOrigin == obj.vectorOrigin) &&
                (minFlag == obj.minFlag) &&
                (maxFlag == obj.maxFlag) &&
                (limitsMode == obj.limitsMode) &&
                (min == obj.min) &&
                (max == obj.max) &&
                (lineStem == obj.lineStem) &&
                (highQuality == obj.highQuality) &&
                (stemWidth == obj.stemWidth));
    }

    public String GetName() { return "Vector"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetUseStride(boolean useStride_)
    {
        useStride = useStride_;
        Select(0);
    }

    public void SetStride(int stride_)
    {
        stride = stride_;
        Select(1);
    }

    public void SetNVectors(int nVectors_)
    {
        nVectors = nVectors_;
        Select(2);
    }

    public void SetLineStyle(int lineStyle_)
    {
        lineStyle = lineStyle_;
        Select(3);
    }

    public void SetLineWidth(int lineWidth_)
    {
        lineWidth = lineWidth_;
        Select(4);
    }

    public void SetScale(double scale_)
    {
        scale = scale_;
        Select(5);
    }

    public void SetScaleByMagnitude(boolean scaleByMagnitude_)
    {
        scaleByMagnitude = scaleByMagnitude_;
        Select(6);
    }

    public void SetAutoScale(boolean autoScale_)
    {
        autoScale = autoScale_;
        Select(7);
    }

    public void SetHeadSize(double headSize_)
    {
        headSize = headSize_;
        Select(8);
    }

    public void SetHeadOn(boolean headOn_)
    {
        headOn = headOn_;
        Select(9);
    }

    public void SetColorByMag(boolean colorByMag_)
    {
        colorByMag = colorByMag_;
        Select(10);
    }

    public void SetUseLegend(boolean useLegend_)
    {
        useLegend = useLegend_;
        Select(11);
    }

    public void SetVectorColor(ColorAttribute vectorColor_)
    {
        vectorColor = vectorColor_;
        Select(12);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(13);
    }

    public void SetVectorOrigin(int vectorOrigin_)
    {
        vectorOrigin = vectorOrigin_;
        Select(14);
    }

    public void SetMinFlag(boolean minFlag_)
    {
        minFlag = minFlag_;
        Select(15);
    }

    public void SetMaxFlag(boolean maxFlag_)
    {
        maxFlag = maxFlag_;
        Select(16);
    }

    public void SetLimitsMode(int limitsMode_)
    {
        limitsMode = limitsMode_;
        Select(17);
    }

    public void SetMin(double min_)
    {
        min = min_;
        Select(18);
    }

    public void SetMax(double max_)
    {
        max = max_;
        Select(19);
    }

    public void SetLineStem(boolean lineStem_)
    {
        lineStem = lineStem_;
        Select(20);
    }

    public void SetHighQuality(boolean highQuality_)
    {
        highQuality = highQuality_;
        Select(21);
    }

    public void SetStemWidth(double stemWidth_)
    {
        stemWidth = stemWidth_;
        Select(22);
    }

    // Property getting methods
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
    public int            GetVectorOrigin() { return vectorOrigin; }
    public boolean        GetMinFlag() { return minFlag; }
    public boolean        GetMaxFlag() { return maxFlag; }
    public int            GetLimitsMode() { return limitsMode; }
    public double         GetMin() { return min; }
    public double         GetMax() { return max; }
    public boolean        GetLineStem() { return lineStem; }
    public boolean        GetHighQuality() { return highQuality; }
    public double         GetStemWidth() { return stemWidth; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(useStride);
        if(WriteSelect(1, buf))
            buf.WriteInt(stride);
        if(WriteSelect(2, buf))
            buf.WriteInt(nVectors);
        if(WriteSelect(3, buf))
            buf.WriteInt(lineStyle);
        if(WriteSelect(4, buf))
            buf.WriteInt(lineWidth);
        if(WriteSelect(5, buf))
            buf.WriteDouble(scale);
        if(WriteSelect(6, buf))
            buf.WriteBool(scaleByMagnitude);
        if(WriteSelect(7, buf))
            buf.WriteBool(autoScale);
        if(WriteSelect(8, buf))
            buf.WriteDouble(headSize);
        if(WriteSelect(9, buf))
            buf.WriteBool(headOn);
        if(WriteSelect(10, buf))
            buf.WriteBool(colorByMag);
        if(WriteSelect(11, buf))
            buf.WriteBool(useLegend);
        if(WriteSelect(12, buf))
            vectorColor.Write(buf);
        if(WriteSelect(13, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(14, buf))
            buf.WriteInt(vectorOrigin);
        if(WriteSelect(15, buf))
            buf.WriteBool(minFlag);
        if(WriteSelect(16, buf))
            buf.WriteBool(maxFlag);
        if(WriteSelect(17, buf))
            buf.WriteInt(limitsMode);
        if(WriteSelect(18, buf))
            buf.WriteDouble(min);
        if(WriteSelect(19, buf))
            buf.WriteDouble(max);
        if(WriteSelect(20, buf))
            buf.WriteBool(lineStem);
        if(WriteSelect(21, buf))
            buf.WriteBool(highQuality);
        if(WriteSelect(22, buf))
            buf.WriteDouble(stemWidth);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetUseStride(buf.ReadBool());
                break;
            case 1:
                SetStride(buf.ReadInt());
                break;
            case 2:
                SetNVectors(buf.ReadInt());
                break;
            case 3:
                SetLineStyle(buf.ReadInt());
                break;
            case 4:
                SetLineWidth(buf.ReadInt());
                break;
            case 5:
                SetScale(buf.ReadDouble());
                break;
            case 6:
                SetScaleByMagnitude(buf.ReadBool());
                break;
            case 7:
                SetAutoScale(buf.ReadBool());
                break;
            case 8:
                SetHeadSize(buf.ReadDouble());
                break;
            case 9:
                SetHeadOn(buf.ReadBool());
                break;
            case 10:
                SetColorByMag(buf.ReadBool());
                break;
            case 11:
                SetUseLegend(buf.ReadBool());
                break;
            case 12:
                vectorColor.Read(buf);
                Select(12);
                break;
            case 13:
                SetColorTableName(buf.ReadString());
                break;
            case 14:
                SetVectorOrigin(buf.ReadInt());
                break;
            case 15:
                SetMinFlag(buf.ReadBool());
                break;
            case 16:
                SetMaxFlag(buf.ReadBool());
                break;
            case 17:
                SetLimitsMode(buf.ReadInt());
                break;
            case 18:
                SetMin(buf.ReadDouble());
                break;
            case 19:
                SetMax(buf.ReadDouble());
                break;
            case 20:
                SetLineStem(buf.ReadBool());
                break;
            case 21:
                SetHighQuality(buf.ReadBool());
                break;
            case 22:
                SetStemWidth(buf.ReadDouble());
                break;
            }
        }
    }


    // Attributes
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
    private int            vectorOrigin;
    private boolean        minFlag;
    private boolean        maxFlag;
    private int            limitsMode;
    private double         min;
    private double         max;
    private boolean        lineStem;
    private boolean        highQuality;
    private double         stemWidth;
}

