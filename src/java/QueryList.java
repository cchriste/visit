package llnl.visit;

import java.util.Vector;
import java.lang.Integer;

// ****************************************************************************
// Class: QueryList
//
// Purpose:
//    List of supported queries
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Tue Dec 2 09:49:27 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

public class QueryList extends AttributeSubject
{
    // Constants
    public final static int QUERYTYPE_DATABASEQUERY = 0;
    public final static int QUERYTYPE_POINTQUERY = 1;
    public final static int QUERYTYPE_LINEQUERY = 2;

    public final static int COORDINATEREPRESENTATION_WORLDSPACE = 0;
    public final static int COORDINATEREPRESENTATION_SCREENSPACE = 1;

    public final static int WINDOWTYPE_BASIC = 0;
    public final static int WINDOWTYPE_SINGLEPOINT = 1;
    public final static int WINDOWTYPE_DOUBLEPOINT = 2;
    public final static int WINDOWTYPE_DOMAINNODE = 3;
    public final static int WINDOWTYPE_DOMAINZONE = 4;
    public final static int WINDOWTYPE_CURRENTPLOT = 5;
    public final static int WINDOWTYPE_CURRENTPLOTVARS = 6;


    public QueryList()
    {
        super(7);

        names = new Vector();
        types = new Vector();
        coordRep = new Vector();
        numInputs = new Vector();
        allowedVarTypes = new Vector();
        winType = new Vector();
        currentPlot = new Vector();
    }

    public QueryList(QueryList obj)
    {
        super(7);

        int i;

        names = new Vector(obj.names.size());
        for(i = 0; i < obj.names.size(); ++i)
            names.addElement(new String((String)obj.names.elementAt(i)));

        types = new Vector();
        for(i = 0; i < obj.types.size(); ++i)
        {
            Integer iv = (Integer)obj.types.elementAt(i);
            types.addElement(new Integer(iv.intValue()));
        }
        coordRep = new Vector();
        for(i = 0; i < obj.coordRep.size(); ++i)
        {
            Integer iv = (Integer)obj.coordRep.elementAt(i);
            coordRep.addElement(new Integer(iv.intValue()));
        }
        numInputs = new Vector();
        for(i = 0; i < obj.numInputs.size(); ++i)
        {
            Integer iv = (Integer)obj.numInputs.elementAt(i);
            numInputs.addElement(new Integer(iv.intValue()));
        }
        allowedVarTypes = new Vector();
        for(i = 0; i < obj.allowedVarTypes.size(); ++i)
        {
            Integer iv = (Integer)obj.allowedVarTypes.elementAt(i);
            allowedVarTypes.addElement(new Integer(iv.intValue()));
        }
        winType = new Vector();
        for(i = 0; i < obj.winType.size(); ++i)
        {
            Integer iv = (Integer)obj.winType.elementAt(i);
            winType.addElement(new Integer(iv.intValue()));
        }
        currentPlot = new Vector();
        for(i = 0; i < obj.currentPlot.size(); ++i)
        {
            Integer iv = (Integer)obj.currentPlot.elementAt(i);
            currentPlot.addElement(new Integer(iv.intValue()));
        }

        SelectAll();
    }

    public boolean equals(QueryList obj)
    {
        int i;

        // Create the return value
        return ((names == obj.names) &&
                (types == obj.types) &&
                (coordRep == obj.coordRep) &&
                (numInputs == obj.numInputs) &&
                (allowedVarTypes == obj.allowedVarTypes) &&
                (winType == obj.winType) &&
                (currentPlot == obj.currentPlot));
    }

    // Property setting methods
    public void SetNames(Vector names_)
    {
        names = names_;
        Select(0);
    }

    public void SetTypes(Vector types_)
    {
        types = types_;
        Select(1);
    }

    public void SetCoordRep(Vector coordRep_)
    {
        coordRep = coordRep_;
        Select(2);
    }

    public void SetNumInputs(Vector numInputs_)
    {
        numInputs = numInputs_;
        Select(3);
    }

    public void SetAllowedVarTypes(Vector allowedVarTypes_)
    {
        allowedVarTypes = allowedVarTypes_;
        Select(4);
    }

    public void SetWinType(Vector winType_)
    {
        winType = winType_;
        Select(5);
    }

    public void SetCurrentPlot(Vector currentPlot_)
    {
        currentPlot = currentPlot_;
        Select(6);
    }

    // Property getting methods
    public Vector GetNames() { return names; }
    public Vector GetTypes() { return types; }
    public Vector GetCoordRep() { return coordRep; }
    public Vector GetNumInputs() { return numInputs; }
    public Vector GetAllowedVarTypes() { return allowedVarTypes; }
    public Vector GetWinType() { return winType; }
    public Vector GetCurrentPlot() { return currentPlot; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteStringVector(names);
        if(WriteSelect(1, buf))
            buf.WriteIntVector(types);
        if(WriteSelect(2, buf))
            buf.WriteIntVector(coordRep);
        if(WriteSelect(3, buf))
            buf.WriteIntVector(numInputs);
        if(WriteSelect(4, buf))
            buf.WriteIntVector(allowedVarTypes);
        if(WriteSelect(5, buf))
            buf.WriteIntVector(winType);
        if(WriteSelect(6, buf))
            buf.WriteIntVector(currentPlot);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetNames(buf.ReadStringVector());
                break;
            case 1:
                SetTypes(buf.ReadIntVector());
                break;
            case 2:
                SetCoordRep(buf.ReadIntVector());
                break;
            case 3:
                SetNumInputs(buf.ReadIntVector());
                break;
            case 4:
                SetAllowedVarTypes(buf.ReadIntVector());
                break;
            case 5:
                SetWinType(buf.ReadIntVector());
                break;
            case 6:
                SetCurrentPlot(buf.ReadIntVector());
                break;
            }
        }
    }


    // Attributes
    private Vector names; // vector of String objects
    private Vector types; // vector of Integer objects
    private Vector coordRep; // vector of Integer objects
    private Vector numInputs; // vector of Integer objects
    private Vector allowedVarTypes; // vector of Integer objects
    private Vector winType; // vector of Integer objects
    private Vector currentPlot; // vector of Integer objects
}

