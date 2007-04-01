package llnl.visit.plots;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;
import llnl.visit.ColorAttribute;

// ****************************************************************************
// Class: TensorAttributes
//
// Purpose:
//    Attributes for the tensor plot
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Nov 24 11:31:17 PDT 2004
//
// Modifications:
//   
// ****************************************************************************

public class TensorAttributes extends AttributeSubject implements Plugin
{
    public TensorAttributes()
    {
        super(10);

        useStride = false;
        stride = 1;
        nTensors = 400;
        scale = 0.25;
        scaleByMagnitude = true;
        autoScale = true;
        colorByEigenvalues = true;
        useLegend = true;
        tensorColor = new ColorAttribute(0, 0, 0);
        colorTableName = new String("Default");
    }

    public TensorAttributes(TensorAttributes obj)
    {
        super(10);

        useStride = obj.useStride;
        stride = obj.stride;
        nTensors = obj.nTensors;
        scale = obj.scale;
        scaleByMagnitude = obj.scaleByMagnitude;
        autoScale = obj.autoScale;
        colorByEigenvalues = obj.colorByEigenvalues;
        useLegend = obj.useLegend;
        tensorColor = new ColorAttribute(obj.tensorColor);
        colorTableName = new String(obj.colorTableName);

        SelectAll();
    }

    public boolean equals(TensorAttributes obj)
    {
        // Create the return value
        return ((useStride == obj.useStride) &&
                (stride == obj.stride) &&
                (nTensors == obj.nTensors) &&
                (scale == obj.scale) &&
                (scaleByMagnitude == obj.scaleByMagnitude) &&
                (autoScale == obj.autoScale) &&
                (colorByEigenvalues == obj.colorByEigenvalues) &&
                (useLegend == obj.useLegend) &&
                (tensorColor == obj.tensorColor) &&
                (colorTableName == obj.colorTableName));
    }

    public String GetName() { return "Tensor"; }
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

    public void SetNTensors(int nTensors_)
    {
        nTensors = nTensors_;
        Select(2);
    }

    public void SetScale(double scale_)
    {
        scale = scale_;
        Select(3);
    }

    public void SetScaleByMagnitude(boolean scaleByMagnitude_)
    {
        scaleByMagnitude = scaleByMagnitude_;
        Select(4);
    }

    public void SetAutoScale(boolean autoScale_)
    {
        autoScale = autoScale_;
        Select(5);
    }

    public void SetColorByEigenvalues(boolean colorByEigenvalues_)
    {
        colorByEigenvalues = colorByEigenvalues_;
        Select(6);
    }

    public void SetUseLegend(boolean useLegend_)
    {
        useLegend = useLegend_;
        Select(7);
    }

    public void SetTensorColor(ColorAttribute tensorColor_)
    {
        tensorColor = tensorColor_;
        Select(8);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(9);
    }

    // Property getting methods
    public boolean        GetUseStride() { return useStride; }
    public int            GetStride() { return stride; }
    public int            GetNTensors() { return nTensors; }
    public double         GetScale() { return scale; }
    public boolean        GetScaleByMagnitude() { return scaleByMagnitude; }
    public boolean        GetAutoScale() { return autoScale; }
    public boolean        GetColorByEigenvalues() { return colorByEigenvalues; }
    public boolean        GetUseLegend() { return useLegend; }
    public ColorAttribute GetTensorColor() { return tensorColor; }
    public String         GetColorTableName() { return colorTableName; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(useStride);
        if(WriteSelect(1, buf))
            buf.WriteInt(stride);
        if(WriteSelect(2, buf))
            buf.WriteInt(nTensors);
        if(WriteSelect(3, buf))
            buf.WriteDouble(scale);
        if(WriteSelect(4, buf))
            buf.WriteBool(scaleByMagnitude);
        if(WriteSelect(5, buf))
            buf.WriteBool(autoScale);
        if(WriteSelect(6, buf))
            buf.WriteBool(colorByEigenvalues);
        if(WriteSelect(7, buf))
            buf.WriteBool(useLegend);
        if(WriteSelect(8, buf))
            tensorColor.Write(buf);
        if(WriteSelect(9, buf))
            buf.WriteString(colorTableName);
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
                SetNTensors(buf.ReadInt());
                break;
            case 3:
                SetScale(buf.ReadDouble());
                break;
            case 4:
                SetScaleByMagnitude(buf.ReadBool());
                break;
            case 5:
                SetAutoScale(buf.ReadBool());
                break;
            case 6:
                SetColorByEigenvalues(buf.ReadBool());
                break;
            case 7:
                SetUseLegend(buf.ReadBool());
                break;
            case 8:
                tensorColor.Read(buf);
                Select(8);
                break;
            case 9:
                SetColorTableName(buf.ReadString());
                break;
            }
        }
    }


    // Attributes
    private boolean        useStride;
    private int            stride;
    private int            nTensors;
    private double         scale;
    private boolean        scaleByMagnitude;
    private boolean        autoScale;
    private boolean        colorByEigenvalues;
    private boolean        useLegend;
    private ColorAttribute tensorColor;
    private String         colorTableName;
}

