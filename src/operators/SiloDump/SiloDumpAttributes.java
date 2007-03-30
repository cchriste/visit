package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: SiloDumpAttributes
//
// Purpose:
//    This class contains attributes for the silo dump operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Aug 15 14:56:28 PST 2002
//
// Modifications:
//   
// ****************************************************************************

public class SiloDumpAttributes extends AttributeSubject implements Plugin
{
    public SiloDumpAttributes()
    {
        super(2);

        filename = new String("dump");
        display = false;
    }

    public SiloDumpAttributes(SiloDumpAttributes obj)
    {
        super(2);

        filename = new String(obj.filename);
        display = obj.display;

        SelectAll();
    }

    public boolean equals(SiloDumpAttributes obj)
    {
        // Create the return value
        return ((filename == obj.filename) &&
                (display == obj.display));
    }

    public String GetName() { return "SiloDump"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetFilename(String filename_)
    {
        filename = filename_;
        Select(0);
    }

    public void SetDisplay(boolean display_)
    {
        display = display_;
        Select(1);
    }

    // Property getting methods
    public String  GetFilename() { return filename; }
    public boolean GetDisplay() { return display; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(filename);
        if(WriteSelect(1, buf))
            buf.WriteBool(display);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetFilename(buf.ReadString());
                break;
            case 1:
                SetDisplay(buf.ReadBool());
                break;
            }
        }
    }


    // Attributes
    private String  filename;
    private boolean display;
}

