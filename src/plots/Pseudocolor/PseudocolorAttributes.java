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
// Creation:   Thu Jul 31 16:07:59 PST 2003
//
// Modifications:
//   
// ****************************************************************************

public class PseudocolorAttributes extends AttributeSubject implements Plugin
{
    // Constants
    public final static int CENTERING_NATURAL = 0;
    public final static int CENTERING_NODAL = 1;
    public final static int CENTERING_ZONAL = 2;

    public final static int SCALING_LINEAR = 0;
    public final static int SCALING_LOG = 1;
    public final static int SCALING_SKEW = 2;

    public final static int LIMITSMODE_ORIGINALDATA = 0;
    public final static int LIMITSMODE_CURRENTPLOT = 1;


    public PseudocolorAttributes()
    {
        super(16);

        legendFlag = true;
        lightingFlag = true;
        minFlag = false;
        maxFlag = false;
        centering = CENTERING_NATURAL;
        scaling = SCALING_LINEAR;
        limitsMode = LIMITSMODE_ORIGINALDATA;
        min = 0;
        max = 1;
        pointSize = 0.05;
        skewFactor = 1;
        opacity = 1;
        colorTableName = new String("hot");
        smoothingLevel = 0;
        pointSizeVarEnabled = false;
        pointSizeVar = new String("default");
    }

    public PseudocolorAttributes(PseudocolorAttributes obj)
    {
        super(16);

        legendFlag = obj.legendFlag;
        lightingFlag = obj.lightingFlag;
        minFlag = obj.minFlag;
        maxFlag = obj.maxFlag;
        centering = obj.centering;
        scaling = obj.scaling;
        limitsMode = obj.limitsMode;
        min = obj.min;
        max = obj.max;
        pointSize = obj.pointSize;
        skewFactor = obj.skewFactor;
        opacity = obj.opacity;
        colorTableName = new String(obj.colorTableName);
        smoothingLevel = obj.smoothingLevel;
        pointSizeVarEnabled = obj.pointSizeVarEnabled;
        pointSizeVar = new String(obj.pointSizeVar);

        SelectAll();
    }

    public boolean equals(PseudocolorAttributes obj)
    {
        // Create the return value
        return ((legendFlag == obj.legendFlag) &&
                (lightingFlag == obj.lightingFlag) &&
                (minFlag == obj.minFlag) &&
                (maxFlag == obj.maxFlag) &&
                (centering == obj.centering) &&
                (scaling == obj.scaling) &&
                (limitsMode == obj.limitsMode) &&
                (min == obj.min) &&
                (max == obj.max) &&
                (pointSize == obj.pointSize) &&
                (skewFactor == obj.skewFactor) &&
                (opacity == obj.opacity) &&
                (colorTableName == obj.colorTableName) &&
                (smoothingLevel == obj.smoothingLevel) &&
                (pointSizeVarEnabled == obj.pointSizeVarEnabled) &&
                (pointSizeVar == obj.pointSizeVar));
    }

    public String GetName() { return "Pseudocolor"; }
    public String GetVersion() { return "1.0"; }

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

    public void SetMinFlag(boolean minFlag_)
    {
        minFlag = minFlag_;
        Select(2);
    }

    public void SetMaxFlag(boolean maxFlag_)
    {
        maxFlag = maxFlag_;
        Select(3);
    }

    public void SetCentering(int centering_)
    {
        centering = centering_;
        Select(4);
    }

    public void SetScaling(int scaling_)
    {
        scaling = scaling_;
        Select(5);
    }

    public void SetLimitsMode(int limitsMode_)
    {
        limitsMode = limitsMode_;
        Select(6);
    }

    public void SetMin(double min_)
    {
        min = min_;
        Select(7);
    }

    public void SetMax(double max_)
    {
        max = max_;
        Select(8);
    }

    public void SetPointSize(double pointSize_)
    {
        pointSize = pointSize_;
        Select(9);
    }

    public void SetSkewFactor(double skewFactor_)
    {
        skewFactor = skewFactor_;
        Select(10);
    }

    public void SetOpacity(double opacity_)
    {
        opacity = opacity_;
        Select(11);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(12);
    }

    public void SetSmoothingLevel(int smoothingLevel_)
    {
        smoothingLevel = smoothingLevel_;
        Select(13);
    }

    public void SetPointSizeVarEnabled(boolean pointSizeVarEnabled_)
    {
        pointSizeVarEnabled = pointSizeVarEnabled_;
        Select(14);
    }

    public void SetPointSizeVar(String pointSizeVar_)
    {
        pointSizeVar = pointSizeVar_;
        Select(15);
    }

    // Property getting methods
    public boolean GetLegendFlag() { return legendFlag; }
    public boolean GetLightingFlag() { return lightingFlag; }
    public boolean GetMinFlag() { return minFlag; }
    public boolean GetMaxFlag() { return maxFlag; }
    public int     GetCentering() { return centering; }
    public int     GetScaling() { return scaling; }
    public int     GetLimitsMode() { return limitsMode; }
    public double  GetMin() { return min; }
    public double  GetMax() { return max; }
    public double  GetPointSize() { return pointSize; }
    public double  GetSkewFactor() { return skewFactor; }
    public double  GetOpacity() { return opacity; }
    public String  GetColorTableName() { return colorTableName; }
    public int     GetSmoothingLevel() { return smoothingLevel; }
    public boolean GetPointSizeVarEnabled() { return pointSizeVarEnabled; }
    public String  GetPointSizeVar() { return pointSizeVar; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(legendFlag);
        if(WriteSelect(1, buf))
            buf.WriteBool(lightingFlag);
        if(WriteSelect(2, buf))
            buf.WriteBool(minFlag);
        if(WriteSelect(3, buf))
            buf.WriteBool(maxFlag);
        if(WriteSelect(4, buf))
            buf.WriteInt(centering);
        if(WriteSelect(5, buf))
            buf.WriteInt(scaling);
        if(WriteSelect(6, buf))
            buf.WriteInt(limitsMode);
        if(WriteSelect(7, buf))
            buf.WriteDouble(min);
        if(WriteSelect(8, buf))
            buf.WriteDouble(max);
        if(WriteSelect(9, buf))
            buf.WriteDouble(pointSize);
        if(WriteSelect(10, buf))
            buf.WriteDouble(skewFactor);
        if(WriteSelect(11, buf))
            buf.WriteDouble(opacity);
        if(WriteSelect(12, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(13, buf))
            buf.WriteInt(smoothingLevel);
        if(WriteSelect(14, buf))
            buf.WriteBool(pointSizeVarEnabled);
        if(WriteSelect(15, buf))
            buf.WriteString(pointSizeVar);
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
                SetMinFlag(buf.ReadBool());
                break;
            case 3:
                SetMaxFlag(buf.ReadBool());
                break;
            case 4:
                SetCentering(buf.ReadInt());
                break;
            case 5:
                SetScaling(buf.ReadInt());
                break;
            case 6:
                SetLimitsMode(buf.ReadInt());
                break;
            case 7:
                SetMin(buf.ReadDouble());
                break;
            case 8:
                SetMax(buf.ReadDouble());
                break;
            case 9:
                SetPointSize(buf.ReadDouble());
                break;
            case 10:
                SetSkewFactor(buf.ReadDouble());
                break;
            case 11:
                SetOpacity(buf.ReadDouble());
                break;
            case 12:
                SetColorTableName(buf.ReadString());
                break;
            case 13:
                SetSmoothingLevel(buf.ReadInt());
                break;
            case 14:
                SetPointSizeVarEnabled(buf.ReadBool());
                break;
            case 15:
                SetPointSizeVar(buf.ReadString());
                break;
            }
        }
    }


    // Attributes
    private boolean legendFlag;
    private boolean lightingFlag;
    private boolean minFlag;
    private boolean maxFlag;
    private int     centering;
    private int     scaling;
    private int     limitsMode;
    private double  min;
    private double  max;
    private double  pointSize;
    private double  skewFactor;
    private double  opacity;
    private String  colorTableName;
    private int     smoothingLevel;
    private boolean pointSizeVarEnabled;
    private String  pointSizeVar;
}

