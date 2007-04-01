package llnl.visit;


// ****************************************************************************
// Class: GlobalLineoutAttributes
//
// Purpose:
//    This file contains global attributes controlling Lineouts.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Feb 3 15:57:15 PST 2005
//
// Modifications:
//   
// ****************************************************************************

public class GlobalLineoutAttributes extends AttributeSubject
{
    // Enum values
    public final static int CURVEOPTIONS_UPDATECURVE = 0;
    public final static int CURVEOPTIONS_CREATECURVE = 1;

    public final static int COLOROPTIONS_REPEATCOLOR = 0;
    public final static int COLOROPTIONS_CREATECOLOR = 1;


    public GlobalLineoutAttributes()
    {
        super(8);

        Dynamic = false;
        createWindow = true;
        windowId = 2;
        samplingOn = false;
        numSamples = 50;
        createReflineLabels = false;
        curveOption = CURVEOPTIONS_UPDATECURVE;
        colorOption = COLOROPTIONS_REPEATCOLOR;
    }

    public GlobalLineoutAttributes(GlobalLineoutAttributes obj)
    {
        super(8);

        Dynamic = obj.Dynamic;
        createWindow = obj.createWindow;
        windowId = obj.windowId;
        samplingOn = obj.samplingOn;
        numSamples = obj.numSamples;
        createReflineLabels = obj.createReflineLabels;
        curveOption = obj.curveOption;
        colorOption = obj.colorOption;

        SelectAll();
    }

    public boolean equals(GlobalLineoutAttributes obj)
    {
        // Create the return value
        return ((Dynamic == obj.Dynamic) &&
                (createWindow == obj.createWindow) &&
                (windowId == obj.windowId) &&
                (samplingOn == obj.samplingOn) &&
                (numSamples == obj.numSamples) &&
                (createReflineLabels == obj.createReflineLabels) &&
                (curveOption == obj.curveOption) &&
                (colorOption == obj.colorOption));
    }

    // Property setting methods
    public void SetDynamic(boolean Dynamic_)
    {
        Dynamic = Dynamic_;
        Select(0);
    }

    public void SetCreateWindow(boolean createWindow_)
    {
        createWindow = createWindow_;
        Select(1);
    }

    public void SetWindowId(int windowId_)
    {
        windowId = windowId_;
        Select(2);
    }

    public void SetSamplingOn(boolean samplingOn_)
    {
        samplingOn = samplingOn_;
        Select(3);
    }

    public void SetNumSamples(int numSamples_)
    {
        numSamples = numSamples_;
        Select(4);
    }

    public void SetCreateReflineLabels(boolean createReflineLabels_)
    {
        createReflineLabels = createReflineLabels_;
        Select(5);
    }

    public void SetCurveOption(int curveOption_)
    {
        curveOption = curveOption_;
        Select(6);
    }

    public void SetColorOption(int colorOption_)
    {
        colorOption = colorOption_;
        Select(7);
    }

    // Property getting methods
    public boolean GetDynamic() { return Dynamic; }
    public boolean GetCreateWindow() { return createWindow; }
    public int     GetWindowId() { return windowId; }
    public boolean GetSamplingOn() { return samplingOn; }
    public int     GetNumSamples() { return numSamples; }
    public boolean GetCreateReflineLabels() { return createReflineLabels; }
    public int     GetCurveOption() { return curveOption; }
    public int     GetColorOption() { return colorOption; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(Dynamic);
        if(WriteSelect(1, buf))
            buf.WriteBool(createWindow);
        if(WriteSelect(2, buf))
            buf.WriteInt(windowId);
        if(WriteSelect(3, buf))
            buf.WriteBool(samplingOn);
        if(WriteSelect(4, buf))
            buf.WriteInt(numSamples);
        if(WriteSelect(5, buf))
            buf.WriteBool(createReflineLabels);
        if(WriteSelect(6, buf))
            buf.WriteInt(curveOption);
        if(WriteSelect(7, buf))
            buf.WriteInt(colorOption);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetDynamic(buf.ReadBool());
                break;
            case 1:
                SetCreateWindow(buf.ReadBool());
                break;
            case 2:
                SetWindowId(buf.ReadInt());
                break;
            case 3:
                SetSamplingOn(buf.ReadBool());
                break;
            case 4:
                SetNumSamples(buf.ReadInt());
                break;
            case 5:
                SetCreateReflineLabels(buf.ReadBool());
                break;
            case 6:
                SetCurveOption(buf.ReadInt());
                break;
            case 7:
                SetColorOption(buf.ReadInt());
                break;
            }
        }
    }


    // Attributes
    private boolean Dynamic;
    private boolean createWindow;
    private int     windowId;
    private boolean samplingOn;
    private int     numSamples;
    private boolean createReflineLabels;
    private int     curveOption;
    private int     colorOption;
}

