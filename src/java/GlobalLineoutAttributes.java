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
// Creation:   Thu Mar 20 10:27:03 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

public class GlobalLineoutAttributes extends AttributeSubject
{
    public GlobalLineoutAttributes()
    {
        super(1);

        Dynamic = false;
    }

    public GlobalLineoutAttributes(GlobalLineoutAttributes obj)
    {
        super(1);

        Dynamic = obj.Dynamic;

        SelectAll();
    }

    public boolean equals(GlobalLineoutAttributes obj)
    {
        // Create the return value
        return ((Dynamic == obj.Dynamic));
    }

    // Property setting methods
    public void SetDynamic(boolean Dynamic_)
    {
        Dynamic = Dynamic_;
        Select(0);
    }

    // Property getting methods
    public boolean GetDynamic() { return Dynamic; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(Dynamic);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        buf.ReadByte();
        SetDynamic(buf.ReadBool());
    }


    // Attributes
    private boolean Dynamic;
}

