package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: InverseGhostZoneAttributes
//
// Purpose:
//    This class contains attributes for the inverse ghost zone operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Jan 8 09:57:13 PDT 2004
//
// Modifications:
//   
// ****************************************************************************

public class InverseGhostZoneAttributes extends AttributeSubject implements Plugin
{
    // Constants
    public final static int SHOWTYPE_GHOSTZONESONLY = 0;
    public final static int SHOWTYPE_GHOSTZONESANDREALZONES = 1;


    public InverseGhostZoneAttributes()
    {
        super(1);

        showType = SHOWTYPE_GHOSTZONESONLY;
    }

    public InverseGhostZoneAttributes(InverseGhostZoneAttributes obj)
    {
        super(1);

        showType = obj.showType;

        SelectAll();
    }

    public boolean equals(InverseGhostZoneAttributes obj)
    {
        // Create the return value
        return ((showType == obj.showType));
    }

    public String GetName() { return "InverseGhostZone"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetShowType(int showType_)
    {
        showType = showType_;
        Select(0);
    }

    // Property getting methods
    public int GetShowType() { return showType; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(showType);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        buf.ReadByte();
        SetShowType(buf.ReadInt());
    }


    // Attributes
    private int showType;
}

