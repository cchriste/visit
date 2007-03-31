package llnl.visit;

import java.util.Vector;

// ****************************************************************************
// Class: GaussianControlPointList
//
// Purpose:
//    This class contains a list of GaussianControlPoint objects.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Jul 31 16:19:15 PST 2003
//
// Modifications:
//   
// ****************************************************************************

public class GaussianControlPointList extends AttributeSubject
{
    public GaussianControlPointList()
    {
        super(1);

        controlPoints = new Vector();
    }

    public GaussianControlPointList(GaussianControlPointList obj)
    {
        super(1);

        int i;

        // *** Copy the controlPoints field ***
        controlPoints = new Vector(obj.controlPoints.size());
        for(i = 0; i < obj.controlPoints.size(); ++i)
        {
            GaussianControlPoint newObj = (GaussianControlPoint)controlPoints.elementAt(i);
            controlPoints.addElement(new GaussianControlPoint(newObj));
        }


        SelectAll();
    }

    public boolean equals(GaussianControlPointList obj)
    {
        int i;

        boolean controlPoints_equal = (obj.controlPoints.size() == controlPoints.size());
        for(i = 0; (i < controlPoints.size()) && controlPoints_equal; ++i)
        {
            // Make references to GaussianControlPoint from Object.
            GaussianControlPoint controlPoints1 = (GaussianControlPoint)controlPoints.elementAt(i);
            GaussianControlPoint controlPoints2 = (GaussianControlPoint)obj.controlPoints.elementAt(i);
            controlPoints_equal = controlPoints1.equals(controlPoints2);
        }

        // Create the return value
        return (controlPoints_equal);
    }

    // Property setting methods
    // Property getting methods
    public Vector GetControlPoints() { return controlPoints; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
        {
            buf.WriteInt(controlPoints.size());
            for(int i = 0; i < controlPoints.size(); ++i)
            {
                GaussianControlPoint tmp = (GaussianControlPoint)controlPoints.elementAt(i);
                tmp.Write(buf);
            }
        }
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        buf.ReadByte();
        {
            int len = buf.ReadInt();
            controlPoints.clear();
            for(int j = 0; j < len; ++j)
            {
                GaussianControlPoint tmp = new GaussianControlPoint();
                tmp.Read(buf);
                controlPoints.addElement(tmp);
            }
        }
        Select(0);
    }

    // Attributegroup convenience methods
    public void AddGaussianControlPoint(GaussianControlPoint obj)
    {
        controlPoints.addElement(new GaussianControlPoint(obj));
        Select(0);
    }

    public void ClearGaussianControlPoints()
    {
        controlPoints.clear();
        Select(0);
    }

    public void RemoveGaussianControlPoint(int index)
    {
        if(index >= 0 && index < controlPoints.size())
        {
            controlPoints.remove(index);
            Select(0);
        }
    }

    public int GetNumGaussianControlPoints()
    {
        return controlPoints.size();
    }

    public GaussianControlPoint GetGaussianControlPoint(int i)
    {
        GaussianControlPoint tmp = (GaussianControlPoint)controlPoints.elementAt(i);
        return tmp;
    }


    // Attributes
    private Vector controlPoints; // vector of GaussianControlPoint objects
}

