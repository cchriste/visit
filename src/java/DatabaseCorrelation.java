package llnl.visit;

import java.util.Vector;
import java.lang.Integer;
import java.lang.Double;

// ****************************************************************************
// Class: DatabaseCorrelation
//
// Purpose:
//    This class encapsulates a database correlation, which is a mapping of one or more databases to a set of indices that go from 0 to N.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Fri Jan 23 14:13:29 PST 2004
//
// Modifications:
//   
// ****************************************************************************

public class DatabaseCorrelation extends AttributeSubject
{
    // Constants
    public final static int CORRELATIONMETHOD_INDEXFORINDEXCORRELATION = 0;
    public final static int CORRELATIONMETHOD_STRETCHEDINDEXCORRELATION = 1;
    public final static int CORRELATIONMETHOD_TIMECORRELATION = 2;
    public final static int CORRELATIONMETHOD_CYCLECORRELATION = 3;
    public final static int CORRELATIONMETHOD_USERDEFINEDCORRELATION = 4;


    public DatabaseCorrelation()
    {
        super(8);

        name = new String("(null)");
        numStates = 1;
        method = CORRELATIONMETHOD_INDEXFORINDEXCORRELATION;
        databaseNames = new Vector();
        databaseNStates = new Vector();
        databaseTimes = new Vector();
        databaseCycles = new Vector();
        indices = new Vector();
    }

    public DatabaseCorrelation(DatabaseCorrelation obj)
    {
        super(8);

        int i;

        name = new String(obj.name);
        numStates = obj.numStates;
        method = obj.method;
        databaseNames = new Vector(obj.databaseNames.size());
        for(i = 0; i < obj.databaseNames.size(); ++i)
            databaseNames.addElement(new String((String)obj.databaseNames.elementAt(i)));

        databaseNStates = new Vector();
        for(i = 0; i < obj.databaseNStates.size(); ++i)
        {
            Integer iv = (Integer)obj.databaseNStates.elementAt(i);
            databaseNStates.addElement(new Integer(iv.intValue()));
        }
        databaseTimes = new Vector(obj.databaseTimes.size());
        for(i = 0; i < obj.databaseTimes.size(); ++i)
        {
            Double dv = (Double)obj.databaseTimes.elementAt(i);
            databaseTimes.addElement(new Double(dv.doubleValue()));
        }

        databaseCycles = new Vector();
        for(i = 0; i < obj.databaseCycles.size(); ++i)
        {
            Integer iv = (Integer)obj.databaseCycles.elementAt(i);
            databaseCycles.addElement(new Integer(iv.intValue()));
        }
        indices = new Vector();
        for(i = 0; i < obj.indices.size(); ++i)
        {
            Integer iv = (Integer)obj.indices.elementAt(i);
            indices.addElement(new Integer(iv.intValue()));
        }

        SelectAll();
    }

    public boolean equals(DatabaseCorrelation obj)
    {
        int i;

        // Create the return value
        return ((name == obj.name) &&
                (numStates == obj.numStates) &&
                (method == obj.method) &&
                (databaseNames == obj.databaseNames) &&
                (databaseNStates == obj.databaseNStates) &&
                (databaseTimes == obj.databaseTimes) &&
                (databaseCycles == obj.databaseCycles) &&
                (indices == obj.indices));
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetNumStates(int numStates_)
    {
        numStates = numStates_;
        Select(1);
    }

    public void SetMethod(int method_)
    {
        method = method_;
        Select(2);
    }

    public void SetDatabaseNames(Vector databaseNames_)
    {
        databaseNames = databaseNames_;
        Select(3);
    }

    public void SetDatabaseNStates(Vector databaseNStates_)
    {
        databaseNStates = databaseNStates_;
        Select(4);
    }

    public void SetDatabaseTimes(Vector databaseTimes_)
    {
        databaseTimes = databaseTimes_;
        Select(5);
    }

    public void SetDatabaseCycles(Vector databaseCycles_)
    {
        databaseCycles = databaseCycles_;
        Select(6);
    }

    public void SetIndices(Vector indices_)
    {
        indices = indices_;
        Select(7);
    }

    // Property getting methods
    public String GetName() { return name; }
    public int    GetNumStates() { return numStates; }
    public int    GetMethod() { return method; }
    public Vector GetDatabaseNames() { return databaseNames; }
    public Vector GetDatabaseNStates() { return databaseNStates; }
    public Vector GetDatabaseTimes() { return databaseTimes; }
    public Vector GetDatabaseCycles() { return databaseCycles; }
    public Vector GetIndices() { return indices; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteInt(numStates);
        if(WriteSelect(2, buf))
            buf.WriteInt(method);
        if(WriteSelect(3, buf))
            buf.WriteStringVector(databaseNames);
        if(WriteSelect(4, buf))
            buf.WriteIntVector(databaseNStates);
        if(WriteSelect(5, buf))
            buf.WriteDoubleVector(databaseTimes);
        if(WriteSelect(6, buf))
            buf.WriteIntVector(databaseCycles);
        if(WriteSelect(7, buf))
            buf.WriteIntVector(indices);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetName(buf.ReadString());
                break;
            case 1:
                SetNumStates(buf.ReadInt());
                break;
            case 2:
                SetMethod(buf.ReadInt());
                break;
            case 3:
                SetDatabaseNames(buf.ReadStringVector());
                break;
            case 4:
                SetDatabaseNStates(buf.ReadIntVector());
                break;
            case 5:
                SetDatabaseTimes(buf.ReadDoubleVector());
                break;
            case 6:
                SetDatabaseCycles(buf.ReadIntVector());
                break;
            case 7:
                SetIndices(buf.ReadIntVector());
                break;
            }
        }
    }


    // Attributes
    private String name;
    private int    numStates;
    private int    method;
    private Vector databaseNames; // vector of String objects
    private Vector databaseNStates; // vector of Integer objects
    private Vector databaseTimes; // vector of Double objects
    private Vector databaseCycles; // vector of Integer objects
    private Vector indices; // vector of Integer objects
}

