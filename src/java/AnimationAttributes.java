package llnl.visit;


// ****************************************************************************
// Class: AnimationAttributes
//
// Purpose:
//    This class contains the animation attributes.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Jul 31 16:19:10 PST 2003
//
// Modifications:
//   
// ****************************************************************************

public class AnimationAttributes extends AttributeSubject
{
    public AnimationAttributes()
    {
        super(2);

        pipelineCachingMode = false;
        timeout = 1;
    }

    public AnimationAttributes(AnimationAttributes obj)
    {
        super(2);

        pipelineCachingMode = obj.pipelineCachingMode;
        timeout = obj.timeout;

        SelectAll();
    }

    public boolean equals(AnimationAttributes obj)
    {
        // Create the return value
        return ((pipelineCachingMode == obj.pipelineCachingMode) &&
                (timeout == obj.timeout));
    }

    // Property setting methods
    public void SetPipelineCachingMode(boolean pipelineCachingMode_)
    {
        pipelineCachingMode = pipelineCachingMode_;
        Select(0);
    }

    public void SetTimeout(int timeout_)
    {
        timeout = timeout_;
        Select(1);
    }

    // Property getting methods
    public boolean GetPipelineCachingMode() { return pipelineCachingMode; }
    public int     GetTimeout() { return timeout; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(pipelineCachingMode);
        if(WriteSelect(1, buf))
            buf.WriteInt(timeout);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetPipelineCachingMode(buf.ReadBool());
                break;
            case 1:
                SetTimeout(buf.ReadInt());
                break;
            }
        }
    }


    // Attributes
    private boolean pipelineCachingMode;
    private int     timeout;
}

