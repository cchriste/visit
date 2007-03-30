package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: BoxAttributes
//
// Purpose:
//    This class contains attributes for the box operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Aug 15 14:56:27 PST 2002
//
// Modifications:
//   
// ****************************************************************************

public class BoxAttributes extends AttributeSubject implements Plugin
{
    // Constants
    public final static int AMOUNT_SOME = 0;
    public final static int AMOUNT_ALL = 1;


    public BoxAttributes()
    {
        super(7);

        amount = AMOUNT_SOME;
        minx = 0;
        maxx = 1;
        miny = 0;
        maxy = 1;
        minz = 0;
        maxz = 1;
    }

    public BoxAttributes(BoxAttributes obj)
    {
        super(7);

        amount = obj.amount;
        minx = obj.minx;
        maxx = obj.maxx;
        miny = obj.miny;
        maxy = obj.maxy;
        minz = obj.minz;
        maxz = obj.maxz;

        SelectAll();
    }

    public boolean equals(BoxAttributes obj)
    {
        // Create the return value
        return ((amount == obj.amount) &&
                (minx == obj.minx) &&
                (maxx == obj.maxx) &&
                (miny == obj.miny) &&
                (maxy == obj.maxy) &&
                (minz == obj.minz) &&
                (maxz == obj.maxz));
    }

    public String GetName() { return "Box"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetAmount(int amount_)
    {
        amount = amount_;
        Select(0);
    }

    public void SetMinx(double minx_)
    {
        minx = minx_;
        Select(1);
    }

    public void SetMaxx(double maxx_)
    {
        maxx = maxx_;
        Select(2);
    }

    public void SetMiny(double miny_)
    {
        miny = miny_;
        Select(3);
    }

    public void SetMaxy(double maxy_)
    {
        maxy = maxy_;
        Select(4);
    }

    public void SetMinz(double minz_)
    {
        minz = minz_;
        Select(5);
    }

    public void SetMaxz(double maxz_)
    {
        maxz = maxz_;
        Select(6);
    }

    // Property getting methods
    public int    GetAmount() { return amount; }
    public double GetMinx() { return minx; }
    public double GetMaxx() { return maxx; }
    public double GetMiny() { return miny; }
    public double GetMaxy() { return maxy; }
    public double GetMinz() { return minz; }
    public double GetMaxz() { return maxz; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(amount);
        if(WriteSelect(1, buf))
            buf.WriteDouble(minx);
        if(WriteSelect(2, buf))
            buf.WriteDouble(maxx);
        if(WriteSelect(3, buf))
            buf.WriteDouble(miny);
        if(WriteSelect(4, buf))
            buf.WriteDouble(maxy);
        if(WriteSelect(5, buf))
            buf.WriteDouble(minz);
        if(WriteSelect(6, buf))
            buf.WriteDouble(maxz);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetAmount(buf.ReadInt());
                break;
            case 1:
                SetMinx(buf.ReadDouble());
                break;
            case 2:
                SetMaxx(buf.ReadDouble());
                break;
            case 3:
                SetMiny(buf.ReadDouble());
                break;
            case 4:
                SetMaxy(buf.ReadDouble());
                break;
            case 5:
                SetMinz(buf.ReadDouble());
                break;
            case 6:
                SetMaxz(buf.ReadDouble());
                break;
            }
        }
    }


    // Attributes
    private int    amount;
    private double minx;
    private double maxx;
    private double miny;
    private double maxy;
    private double minz;
    private double maxz;
}

