package llnl.visit;

import java.util.Vector;
import java.lang.Integer;

// ****************************************************************************
// Class: EngineList
//
// Purpose:
//    This class contains a list of host names on which engines are running.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Tue Mar 30 17:34:14 PST 2004
//
// Modifications:
//   
// ****************************************************************************

public class EngineList extends AttributeSubject
{
    public EngineList()
    {
        super(5);

        engines = new Vector();
        numProcessors = new Vector();
        numNodes = new Vector();
        loadBalancing = new Vector();
        simulationName = new Vector();
    }

    public EngineList(EngineList obj)
    {
        super(5);

        int i;

        engines = new Vector(obj.engines.size());
        for(i = 0; i < obj.engines.size(); ++i)
            engines.addElement(new String((String)obj.engines.elementAt(i)));

        numProcessors = new Vector();
        for(i = 0; i < obj.numProcessors.size(); ++i)
        {
            Integer iv = (Integer)obj.numProcessors.elementAt(i);
            numProcessors.addElement(new Integer(iv.intValue()));
        }
        numNodes = new Vector();
        for(i = 0; i < obj.numNodes.size(); ++i)
        {
            Integer iv = (Integer)obj.numNodes.elementAt(i);
            numNodes.addElement(new Integer(iv.intValue()));
        }
        loadBalancing = new Vector();
        for(i = 0; i < obj.loadBalancing.size(); ++i)
        {
            Integer iv = (Integer)obj.loadBalancing.elementAt(i);
            loadBalancing.addElement(new Integer(iv.intValue()));
        }
        simulationName = new Vector(obj.simulationName.size());
        for(i = 0; i < obj.simulationName.size(); ++i)
            simulationName.addElement(new String((String)obj.simulationName.elementAt(i)));


        SelectAll();
    }

    public boolean equals(EngineList obj)
    {
        int i;

        // Create the return value
        return ((engines == obj.engines) &&
                (numProcessors == obj.numProcessors) &&
                (numNodes == obj.numNodes) &&
                (loadBalancing == obj.loadBalancing) &&
                (simulationName == obj.simulationName));
    }

    // Property setting methods
    public void SetEngines(Vector engines_)
    {
        engines = engines_;
        Select(0);
    }

    public void SetNumProcessors(Vector numProcessors_)
    {
        numProcessors = numProcessors_;
        Select(1);
    }

    public void SetNumNodes(Vector numNodes_)
    {
        numNodes = numNodes_;
        Select(2);
    }

    public void SetLoadBalancing(Vector loadBalancing_)
    {
        loadBalancing = loadBalancing_;
        Select(3);
    }

    public void SetSimulationName(Vector simulationName_)
    {
        simulationName = simulationName_;
        Select(4);
    }

    // Property getting methods
    public Vector GetEngines() { return engines; }
    public Vector GetNumProcessors() { return numProcessors; }
    public Vector GetNumNodes() { return numNodes; }
    public Vector GetLoadBalancing() { return loadBalancing; }
    public Vector GetSimulationName() { return simulationName; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteStringVector(engines);
        if(WriteSelect(1, buf))
            buf.WriteIntVector(numProcessors);
        if(WriteSelect(2, buf))
            buf.WriteIntVector(numNodes);
        if(WriteSelect(3, buf))
            buf.WriteIntVector(loadBalancing);
        if(WriteSelect(4, buf))
            buf.WriteStringVector(simulationName);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetEngines(buf.ReadStringVector());
                break;
            case 1:
                SetNumProcessors(buf.ReadIntVector());
                break;
            case 2:
                SetNumNodes(buf.ReadIntVector());
                break;
            case 3:
                SetLoadBalancing(buf.ReadIntVector());
                break;
            case 4:
                SetSimulationName(buf.ReadStringVector());
                break;
            }
        }
    }


    // Attributes
    private Vector engines; // vector of String objects
    private Vector numProcessors; // vector of Integer objects
    private Vector numNodes; // vector of Integer objects
    private Vector loadBalancing; // vector of Integer objects
    private Vector simulationName; // vector of String objects
}

