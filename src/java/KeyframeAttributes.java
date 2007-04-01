package llnl.visit;


// ****************************************************************************
// Class: KeyframeAttributes
//
// Purpose:
//    This class contains the attributes used for keyframing.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Tue Apr 6 23:40:05 PST 2004
//
// Modifications:
//   
// ****************************************************************************

public class KeyframeAttributes extends AttributeSubject
{
    public KeyframeAttributes()
    {
        super(3);

        enabled = false;
        nFrames = 0;
        nFramesWasUserSet = false;
    }

    public KeyframeAttributes(KeyframeAttributes obj)
    {
        super(3);

        enabled = obj.enabled;
        nFrames = obj.nFrames;
        nFramesWasUserSet = obj.nFramesWasUserSet;

        SelectAll();
    }

    public boolean equals(KeyframeAttributes obj)
    {
        // Create the return value
        return ((enabled == obj.enabled) &&
                (nFrames == obj.nFrames) &&
                (nFramesWasUserSet == obj.nFramesWasUserSet));
    }

    // Property setting methods
    public void SetEnabled(boolean enabled_)
    {
        enabled = enabled_;
        Select(0);
    }

    public void SetNFrames(int nFrames_)
    {
        nFrames = nFrames_;
        Select(1);
    }

    public void SetNFramesWasUserSet(boolean nFramesWasUserSet_)
    {
        nFramesWasUserSet = nFramesWasUserSet_;
        Select(2);
    }

    // Property getting methods
    public boolean GetEnabled() { return enabled; }
    public int     GetNFrames() { return nFrames; }
    public boolean GetNFramesWasUserSet() { return nFramesWasUserSet; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(enabled);
        if(WriteSelect(1, buf))
            buf.WriteInt(nFrames);
        if(WriteSelect(2, buf))
            buf.WriteBool(nFramesWasUserSet);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetEnabled(buf.ReadBool());
                break;
            case 1:
                SetNFrames(buf.ReadInt());
                break;
            case 2:
                SetNFramesWasUserSet(buf.ReadBool());
                break;
            }
        }
    }


    // Attributes
    private boolean enabled;
    private int     nFrames;
    private boolean nFramesWasUserSet;
}

