package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: IndexSelectAttributes
//
// Purpose:
//    This class contains attributes for the index select operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Aug 15 14:56:27 PST 2002
//
// Modifications:
//   
// ****************************************************************************

public class IndexSelectAttributes extends AttributeSubject implements Plugin
{
    // Constants
    public final static int DIMENSION_ONED = 0;
    public final static int DIMENSION_TWOD = 1;
    public final static int DIMENSION_THREED = 2;

    public final static int DATATYPE_ALLDOMAINS = 0;
    public final static int DATATYPE_ONEDOMAIN = 1;
    public final static int DATATYPE_ONEGROUP = 2;


    public IndexSelectAttributes()
    {
        super(13);

        dim = DIMENSION_TWOD;
        xMin = 0;
        xMax = -1;
        xIncr = 1;
        yMin = 0;
        yMax = -1;
        yIncr = 1;
        zMin = 0;
        zMax = -1;
        zIncr = 1;
        whichData = DATATYPE_ALLDOMAINS;
        domainIndex = 0;
        groupIndex = 0;
    }

    public IndexSelectAttributes(IndexSelectAttributes obj)
    {
        super(13);

        dim = obj.dim;
        xMin = obj.xMin;
        xMax = obj.xMax;
        xIncr = obj.xIncr;
        yMin = obj.yMin;
        yMax = obj.yMax;
        yIncr = obj.yIncr;
        zMin = obj.zMin;
        zMax = obj.zMax;
        zIncr = obj.zIncr;
        whichData = obj.whichData;
        domainIndex = obj.domainIndex;
        groupIndex = obj.groupIndex;

        SelectAll();
    }

    public boolean equals(IndexSelectAttributes obj)
    {
        // Create the return value
        return ((dim == obj.dim) &&
                (xMin == obj.xMin) &&
                (xMax == obj.xMax) &&
                (xIncr == obj.xIncr) &&
                (yMin == obj.yMin) &&
                (yMax == obj.yMax) &&
                (yIncr == obj.yIncr) &&
                (zMin == obj.zMin) &&
                (zMax == obj.zMax) &&
                (zIncr == obj.zIncr) &&
                (whichData == obj.whichData) &&
                (domainIndex == obj.domainIndex) &&
                (groupIndex == obj.groupIndex));
    }

    public String GetName() { return "IndexSelect"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetDim(int dim_)
    {
        dim = dim_;
        Select(0);
    }

    public void SetXMin(int xMin_)
    {
        xMin = xMin_;
        Select(1);
    }

    public void SetXMax(int xMax_)
    {
        xMax = xMax_;
        Select(2);
    }

    public void SetXIncr(int xIncr_)
    {
        xIncr = xIncr_;
        Select(3);
    }

    public void SetYMin(int yMin_)
    {
        yMin = yMin_;
        Select(4);
    }

    public void SetYMax(int yMax_)
    {
        yMax = yMax_;
        Select(5);
    }

    public void SetYIncr(int yIncr_)
    {
        yIncr = yIncr_;
        Select(6);
    }

    public void SetZMin(int zMin_)
    {
        zMin = zMin_;
        Select(7);
    }

    public void SetZMax(int zMax_)
    {
        zMax = zMax_;
        Select(8);
    }

    public void SetZIncr(int zIncr_)
    {
        zIncr = zIncr_;
        Select(9);
    }

    public void SetWhichData(int whichData_)
    {
        whichData = whichData_;
        Select(10);
    }

    public void SetDomainIndex(int domainIndex_)
    {
        domainIndex = domainIndex_;
        Select(11);
    }

    public void SetGroupIndex(int groupIndex_)
    {
        groupIndex = groupIndex_;
        Select(12);
    }

    // Property getting methods
    public int GetDim() { return dim; }
    public int GetXMin() { return xMin; }
    public int GetXMax() { return xMax; }
    public int GetXIncr() { return xIncr; }
    public int GetYMin() { return yMin; }
    public int GetYMax() { return yMax; }
    public int GetYIncr() { return yIncr; }
    public int GetZMin() { return zMin; }
    public int GetZMax() { return zMax; }
    public int GetZIncr() { return zIncr; }
    public int GetWhichData() { return whichData; }
    public int GetDomainIndex() { return domainIndex; }
    public int GetGroupIndex() { return groupIndex; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(dim);
        if(WriteSelect(1, buf))
            buf.WriteInt(xMin);
        if(WriteSelect(2, buf))
            buf.WriteInt(xMax);
        if(WriteSelect(3, buf))
            buf.WriteInt(xIncr);
        if(WriteSelect(4, buf))
            buf.WriteInt(yMin);
        if(WriteSelect(5, buf))
            buf.WriteInt(yMax);
        if(WriteSelect(6, buf))
            buf.WriteInt(yIncr);
        if(WriteSelect(7, buf))
            buf.WriteInt(zMin);
        if(WriteSelect(8, buf))
            buf.WriteInt(zMax);
        if(WriteSelect(9, buf))
            buf.WriteInt(zIncr);
        if(WriteSelect(10, buf))
            buf.WriteInt(whichData);
        if(WriteSelect(11, buf))
            buf.WriteInt(domainIndex);
        if(WriteSelect(12, buf))
            buf.WriteInt(groupIndex);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetDim(buf.ReadInt());
                break;
            case 1:
                SetXMin(buf.ReadInt());
                break;
            case 2:
                SetXMax(buf.ReadInt());
                break;
            case 3:
                SetXIncr(buf.ReadInt());
                break;
            case 4:
                SetYMin(buf.ReadInt());
                break;
            case 5:
                SetYMax(buf.ReadInt());
                break;
            case 6:
                SetYIncr(buf.ReadInt());
                break;
            case 7:
                SetZMin(buf.ReadInt());
                break;
            case 8:
                SetZMax(buf.ReadInt());
                break;
            case 9:
                SetZIncr(buf.ReadInt());
                break;
            case 10:
                SetWhichData(buf.ReadInt());
                break;
            case 11:
                SetDomainIndex(buf.ReadInt());
                break;
            case 12:
                SetGroupIndex(buf.ReadInt());
                break;
            }
        }
    }


    // Attributes
    private int dim;
    private int xMin;
    private int xMax;
    private int xIncr;
    private int yMin;
    private int yMax;
    private int yIncr;
    private int zMin;
    private int zMax;
    private int zIncr;
    private int whichData;
    private int domainIndex;
    private int groupIndex;
}

