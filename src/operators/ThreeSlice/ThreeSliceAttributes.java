package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: ThreeSliceAttributes
//
// Purpose:
//    This class contains attributes for the threeslice operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Jul 31 16:09:31 PST 2003
//
// Modifications:
//   
// ****************************************************************************

public class ThreeSliceAttributes extends AttributeSubject implements Plugin
{
    public ThreeSliceAttributes()
    {
        super(4);

        x = 0f;
        y = 0f;
        z = 0f;
        interactive = true;
    }

    public ThreeSliceAttributes(ThreeSliceAttributes obj)
    {
        super(4);

        x = obj.x;
        y = obj.y;
        z = obj.z;
        interactive = obj.interactive;

        SelectAll();
    }

    public boolean equals(ThreeSliceAttributes obj)
    {
        // Create the return value
        return ((x == obj.x) &&
                (y == obj.y) &&
                (z == obj.z) &&
                (interactive == obj.interactive));
    }

    public String GetName() { return "ThreeSlice"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetX(float x_)
    {
        x = x_;
        Select(0);
    }

    public void SetY(float y_)
    {
        y = y_;
        Select(1);
    }

    public void SetZ(float z_)
    {
        z = z_;
        Select(2);
    }

    public void SetInteractive(boolean interactive_)
    {
        interactive = interactive_;
        Select(3);
    }

    // Property getting methods
    public float   GetX() { return x; }
    public float   GetY() { return y; }
    public float   GetZ() { return z; }
    public boolean GetInteractive() { return interactive; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteFloat(x);
        if(WriteSelect(1, buf))
            buf.WriteFloat(y);
        if(WriteSelect(2, buf))
            buf.WriteFloat(z);
        if(WriteSelect(3, buf))
            buf.WriteBool(interactive);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetX(buf.ReadFloat());
                break;
            case 1:
                SetY(buf.ReadFloat());
                break;
            case 2:
                SetZ(buf.ReadFloat());
                break;
            case 3:
                SetInteractive(buf.ReadBool());
                break;
            }
        }
    }


    // Attributes
    private float   x;
    private float   y;
    private float   z;
    private boolean interactive;
}

