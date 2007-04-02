package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;
import java.util.Vector;

// ****************************************************************************
// Class: DeferExpressionAttributes
//
// Purpose:
//    Attributes for the DeferExpression operator
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Tue Sep 20 13:25:16 PST 2005
//
// Modifications:
//   
// ****************************************************************************

public class DeferExpressionAttributes extends AttributeSubject implements Plugin
{
    public DeferExpressionAttributes()
    {
        super(1);

        exprs = new Vector();
    }

    public DeferExpressionAttributes(DeferExpressionAttributes obj)
    {
        super(1);

        int i;

        exprs = new Vector(obj.exprs.size());
        for(i = 0; i < obj.exprs.size(); ++i)
            exprs.addElement(new String((String)obj.exprs.elementAt(i)));


        SelectAll();
    }

    public boolean equals(DeferExpressionAttributes obj)
    {
        int i;

        // Create the return value
        return ((exprs == obj.exprs));
    }

    public String GetName() { return "DeferExpression"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetExprs(Vector exprs_)
    {
        exprs = exprs_;
        Select(0);
    }

    // Property getting methods
    public Vector GetExprs() { return exprs; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteStringVector(exprs);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        buf.ReadByte();
        SetExprs(buf.ReadStringVector());
    }


    // Attributes
    private Vector exprs; // vector of String objects
}

