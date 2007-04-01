package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: LineoutAttributes
//
// Purpose:
//    Attributes for the Lineout operator 
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Jul 28 14:34:41 PST 2004
//
// Modifications:
//   
// ****************************************************************************

public class LineoutAttributes extends AttributeSubject implements Plugin
{
    public LineoutAttributes()
    {
        super(7);

        point1 = new double[3];
        point1[0] = 0;
        point1[1] = 0;
        point1[2] = 0;
        point2 = new double[3];
        point2[0] = 1;
        point2[1] = 1;
        point2[2] = 0;
        samplingOn = false;
        numberOfSamplePoints = 50;
        interactive = false;
        reflineLabels = false;
        designator = new String("(null)");
    }

    public LineoutAttributes(LineoutAttributes obj)
    {
        super(7);

        int i;

        point1 = new double[3];
        point1[0] = obj.point1[0];
        point1[1] = obj.point1[1];
        point1[2] = obj.point1[2];

        point2 = new double[3];
        point2[0] = obj.point2[0];
        point2[1] = obj.point2[1];
        point2[2] = obj.point2[2];

        samplingOn = obj.samplingOn;
        numberOfSamplePoints = obj.numberOfSamplePoints;
        interactive = obj.interactive;
        reflineLabels = obj.reflineLabels;
        designator = new String(obj.designator);

        SelectAll();
    }

    public boolean equals(LineoutAttributes obj)
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
                (samplingOn == obj.samplingOn) &&
                (numberOfSamplePoints == obj.numberOfSamplePoints) &&
                (interactive == obj.interactive) &&
                (reflineLabels == obj.reflineLabels) &&
                (designator == obj.designator));
    }

    public String GetName() { return "Lineout"; }
    public String GetVersion() { return "1.0"; }

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

    public void SetSamplingOn(boolean samplingOn_)
    {
        samplingOn = samplingOn_;
        Select(2);
    }

    public void SetNumberOfSamplePoints(int numberOfSamplePoints_)
    {
        numberOfSamplePoints = numberOfSamplePoints_;
        Select(3);
    }

    public void SetInteractive(boolean interactive_)
    {
        interactive = interactive_;
        Select(4);
    }

    public void SetReflineLabels(boolean reflineLabels_)
    {
        reflineLabels = reflineLabels_;
        Select(5);
    }

    public void SetDesignator(String designator_)
    {
        designator = designator_;
        Select(6);
    }

    // Property getting methods
    public double[] GetPoint1() { return point1; }
    public double[] GetPoint2() { return point2; }
    public boolean  GetSamplingOn() { return samplingOn; }
    public int      GetNumberOfSamplePoints() { return numberOfSamplePoints; }
    public boolean  GetInteractive() { return interactive; }
    public boolean  GetReflineLabels() { return reflineLabels; }
    public String   GetDesignator() { return designator; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(point1);
        if(WriteSelect(1, buf))
            buf.WriteDoubleArray(point2);
        if(WriteSelect(2, buf))
            buf.WriteBool(samplingOn);
        if(WriteSelect(3, buf))
            buf.WriteInt(numberOfSamplePoints);
        if(WriteSelect(4, buf))
            buf.WriteBool(interactive);
        if(WriteSelect(5, buf))
            buf.WriteBool(reflineLabels);
        if(WriteSelect(6, buf))
            buf.WriteString(designator);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetPoint1(buf.ReadDoubleArray());
                break;
            case 1:
                SetPoint2(buf.ReadDoubleArray());
                break;
            case 2:
                SetSamplingOn(buf.ReadBool());
                break;
            case 3:
                SetNumberOfSamplePoints(buf.ReadInt());
                break;
            case 4:
                SetInteractive(buf.ReadBool());
                break;
            case 5:
                SetReflineLabels(buf.ReadBool());
                break;
            case 6:
                SetDesignator(buf.ReadString());
                break;
            }
        }
    }


    // Attributes
    private double[] point1;
    private double[] point2;
    private boolean  samplingOn;
    private int      numberOfSamplePoints;
    private boolean  interactive;
    private boolean  reflineLabels;
    private String   designator;
}

