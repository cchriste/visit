package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: DisplaceAttributes
//
// Purpose:
//    This class contains attributes for the displace operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Oct 6 16:08:36 PST 2005
//
// Modifications:
//   
// ****************************************************************************

public class DisplaceAttributes extends AttributeSubject implements Plugin
{
    public DisplaceAttributes()
    {
        super(2);

        factor = 1;
        variable = new String("default");
    }

    public DisplaceAttributes(DisplaceAttributes obj)
    {
        super(2);

        factor = obj.factor;
        variable = new String(obj.variable);

        SelectAll();
    }

    public boolean equals(DisplaceAttributes obj)
    {
        // Create the return value
        return ((factor == obj.factor) &&
                (variable == obj.variable));
    }

    public String GetName() { return "Displace"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetFactor(double factor_)
    {
        factor = factor_;
        Select(0);
    }

    public void SetVariable(String variable_)
    {
        variable = variable_;
        Select(1);
    }

    // Property getting methods
    public double GetFactor() { return factor; }
    public String GetVariable() { return variable; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDouble(factor);
        if(WriteSelect(1, buf))
            buf.WriteString(variable);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetFactor(buf.ReadDouble());
                break;
            case 1:
                SetVariable(buf.ReadString());
                break;
            }
        }
    }


    // Attributes
    private double factor;
    private String variable;
}

