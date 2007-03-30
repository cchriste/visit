package llnl.visit;

import java.lang.Integer;
import java.util.Vector;

// ****************************************************************************
// Class: SILMatrixAttributes
//
// Purpose:
//    This class contain the information needed to represent a SIL Matrix.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Tue Feb 11 07:24:11 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

public class SILMatrixAttributes extends AttributeSubject
{
    public SILMatrixAttributes()
    {
        super(6);

        set1 = new Vector();
        category1 = new String("(null)");
        role1 = 0;
        set2 = new Vector();
        category2 = new String("(null)");
        role2 = 0;
    }

    public SILMatrixAttributes(SILMatrixAttributes obj)
    {
        super(6);

        int i;

        set1 = new Vector();
        for(i = 0; i < obj.set1.size(); ++i)
        {
            Integer iv = (Integer)obj.set1.elementAt(i);
            set1.addElement(new Integer(iv.intValue()));
        }
        category1 = new String(obj.category1);
        role1 = obj.role1;
        set2 = new Vector();
        for(i = 0; i < obj.set2.size(); ++i)
        {
            Integer iv = (Integer)obj.set2.elementAt(i);
            set2.addElement(new Integer(iv.intValue()));
        }
        category2 = new String(obj.category2);
        role2 = obj.role2;

        SelectAll();
    }

    public boolean equals(SILMatrixAttributes obj)
    {
        int i;

        // Create the return value
        return ((set1 == obj.set1) &&
                (category1 == obj.category1) &&
                (role1 == obj.role1) &&
                (set2 == obj.set2) &&
                (category2 == obj.category2) &&
                (role2 == obj.role2));
    }

    // Property setting methods
    public void SetSet1(Vector set1_)
    {
        set1 = set1_;
        Select(0);
    }

    public void SetCategory1(String category1_)
    {
        category1 = category1_;
        Select(1);
    }

    public void SetRole1(int role1_)
    {
        role1 = role1_;
        Select(2);
    }

    public void SetSet2(Vector set2_)
    {
        set2 = set2_;
        Select(3);
    }

    public void SetCategory2(String category2_)
    {
        category2 = category2_;
        Select(4);
    }

    public void SetRole2(int role2_)
    {
        role2 = role2_;
        Select(5);
    }

    // Property getting methods
    public Vector GetSet1() { return set1; }
    public String GetCategory1() { return category1; }
    public int    GetRole1() { return role1; }
    public Vector GetSet2() { return set2; }
    public String GetCategory2() { return category2; }
    public int    GetRole2() { return role2; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteIntVector(set1);
        if(WriteSelect(1, buf))
            buf.WriteString(category1);
        if(WriteSelect(2, buf))
            buf.WriteInt(role1);
        if(WriteSelect(3, buf))
            buf.WriteIntVector(set2);
        if(WriteSelect(4, buf))
            buf.WriteString(category2);
        if(WriteSelect(5, buf))
            buf.WriteInt(role2);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetSet1(buf.ReadIntVector());
                break;
            case 1:
                SetCategory1(buf.ReadString());
                break;
            case 2:
                SetRole1(buf.ReadInt());
                break;
            case 3:
                SetSet2(buf.ReadIntVector());
                break;
            case 4:
                SetCategory2(buf.ReadString());
                break;
            case 5:
                SetRole2(buf.ReadInt());
                break;
            }
        }
    }


    // Attributes
    private Vector set1; // vector of Integer objects
    private String category1;
    private int    role1;
    private Vector set2; // vector of Integer objects
    private String category2;
    private int    role2;
}

