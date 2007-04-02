package llnl.visit;

import java.lang.Integer;
import java.util.Vector;
import java.lang.Double;

// ****************************************************************************
// Class: ClientMethod
//
// Purpose:
//    This class contains the attributes needed to make one client execute a method on another client.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu May 5 18:50:15 PST 2005
//
// Modifications:
//   
// ****************************************************************************

public class ClientMethod extends AttributeSubject
{
    public ClientMethod()
    {
        super(4);

        methodName = new String("");
        intArgs = new Vector();
        doubleArgs = new Vector();
        stringArgs = new Vector();
    }

    public ClientMethod(ClientMethod obj)
    {
        super(4);

        int i;

        methodName = new String(obj.methodName);
        intArgs = new Vector();
        for(i = 0; i < obj.intArgs.size(); ++i)
        {
            Integer iv = (Integer)obj.intArgs.elementAt(i);
            intArgs.addElement(new Integer(iv.intValue()));
        }
        doubleArgs = new Vector(obj.doubleArgs.size());
        for(i = 0; i < obj.doubleArgs.size(); ++i)
        {
            Double dv = (Double)obj.doubleArgs.elementAt(i);
            doubleArgs.addElement(new Double(dv.doubleValue()));
        }

        stringArgs = new Vector(obj.stringArgs.size());
        for(i = 0; i < obj.stringArgs.size(); ++i)
            stringArgs.addElement(new String((String)obj.stringArgs.elementAt(i)));


        SelectAll();
    }

    public boolean equals(ClientMethod obj)
    {
        int i;

        // Create the return value
        return ((methodName == obj.methodName) &&
                (intArgs == obj.intArgs) &&
                (doubleArgs == obj.doubleArgs) &&
                (stringArgs == obj.stringArgs));
    }

    // Property setting methods
    public void SetMethodName(String methodName_)
    {
        methodName = methodName_;
        Select(0);
    }

    public void SetIntArgs(Vector intArgs_)
    {
        intArgs = intArgs_;
        Select(1);
    }

    public void SetDoubleArgs(Vector doubleArgs_)
    {
        doubleArgs = doubleArgs_;
        Select(2);
    }

    public void SetStringArgs(Vector stringArgs_)
    {
        stringArgs = stringArgs_;
        Select(3);
    }

    // Property getting methods
    public String GetMethodName() { return methodName; }
    public Vector GetIntArgs() { return intArgs; }
    public Vector GetDoubleArgs() { return doubleArgs; }
    public Vector GetStringArgs() { return stringArgs; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(methodName);
        if(WriteSelect(1, buf))
            buf.WriteIntVector(intArgs);
        if(WriteSelect(2, buf))
            buf.WriteDoubleVector(doubleArgs);
        if(WriteSelect(3, buf))
            buf.WriteStringVector(stringArgs);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetMethodName(buf.ReadString());
                break;
            case 1:
                SetIntArgs(buf.ReadIntVector());
                break;
            case 2:
                SetDoubleArgs(buf.ReadDoubleVector());
                break;
            case 3:
                SetStringArgs(buf.ReadStringVector());
                break;
            }
        }
    }


    // Attributes
    private String methodName;
    private Vector intArgs; // vector of Integer objects
    private Vector doubleArgs; // vector of Double objects
    private Vector stringArgs; // vector of String objects
}

