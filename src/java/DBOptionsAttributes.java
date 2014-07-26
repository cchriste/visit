// ***************************************************************************
//
// Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-442911
// All rights reserved.
//
// This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
// full copyright notice is contained in the file COPYRIGHT located at the root
// of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
//
// Redistribution  and  use  in  source  and  binary  forms,  with  or  without
// modification, are permitted provided that the following conditions are met:
//
//  - Redistributions of  source code must  retain the above  copyright notice,
//    this list of conditions and the disclaimer below.
//  - Redistributions in binary form must reproduce the above copyright notice,
//    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
//    documentation and/or other materials provided with the distribution.
//  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
//    be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
// LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
// DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
// SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
// CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
// LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
// OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ***************************************************************************

package llnl.visit;

import java.lang.Integer;
import java.util.Vector;
import java.lang.Double;

// ****************************************************************************
// Class: DBOptionsAttributes
//
// Purpose:
//    Attributes of database options
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class DBOptionsAttributes extends AttributeSubject
{
    private static int DBOptionsAttributes_numAdditionalAtts = 11;

    // Enum values
    public final static int OPTIONTYPE_BOOL = 0;
    public final static int OPTIONTYPE_INT = 1;
    public final static int OPTIONTYPE_FLOAT = 2;
    public final static int OPTIONTYPE_DOUBLE = 3;
    public final static int OPTIONTYPE_STRING = 4;
    public final static int OPTIONTYPE_ENUM = 5;


    public DBOptionsAttributes()
    {
        super(DBOptionsAttributes_numAdditionalAtts);

        types = new Vector();
        names = new Vector();
        optBools = new Vector();
        optFloats = new Vector();
        optDoubles = new Vector();
        optInts = new Vector();
        optStrings = new Vector();
        optEnums = new Vector();
        enumStrings = new Vector();
        enumStringsSizes = new Vector();
        obsoleteNames = new Vector();
    }

    public DBOptionsAttributes(int nMoreFields)
    {
        super(DBOptionsAttributes_numAdditionalAtts + nMoreFields);

        types = new Vector();
        names = new Vector();
        optBools = new Vector();
        optFloats = new Vector();
        optDoubles = new Vector();
        optInts = new Vector();
        optStrings = new Vector();
        optEnums = new Vector();
        enumStrings = new Vector();
        enumStringsSizes = new Vector();
        obsoleteNames = new Vector();
    }

    public DBOptionsAttributes(DBOptionsAttributes obj)
    {
        super(DBOptionsAttributes_numAdditionalAtts);

        int i;

        types = new Vector();
        for(i = 0; i < obj.types.size(); ++i)
        {
            Integer iv = (Integer)obj.types.elementAt(i);
            types.addElement(new Integer(iv.intValue()));
        }
        names = new Vector(obj.names.size());
        for(i = 0; i < obj.names.size(); ++i)
            names.addElement(new String((String)obj.names.elementAt(i)));

        optBools = new Vector();
        for(i = 0; i < obj.optBools.size(); ++i)
        {
            Integer iv = (Integer)obj.optBools.elementAt(i);
            optBools.addElement(new Integer(iv.intValue()));
        }
        optFloats = new Vector(obj.optFloats.size());
        for(i = 0; i < obj.optFloats.size(); ++i)
        {
            Double dv = (Double)obj.optFloats.elementAt(i);
            optFloats.addElement(new Double(dv.doubleValue()));
        }

        optDoubles = new Vector(obj.optDoubles.size());
        for(i = 0; i < obj.optDoubles.size(); ++i)
        {
            Double dv = (Double)obj.optDoubles.elementAt(i);
            optDoubles.addElement(new Double(dv.doubleValue()));
        }

        optInts = new Vector();
        for(i = 0; i < obj.optInts.size(); ++i)
        {
            Integer iv = (Integer)obj.optInts.elementAt(i);
            optInts.addElement(new Integer(iv.intValue()));
        }
        optStrings = new Vector(obj.optStrings.size());
        for(i = 0; i < obj.optStrings.size(); ++i)
            optStrings.addElement(new String((String)obj.optStrings.elementAt(i)));

        optEnums = new Vector();
        for(i = 0; i < obj.optEnums.size(); ++i)
        {
            Integer iv = (Integer)obj.optEnums.elementAt(i);
            optEnums.addElement(new Integer(iv.intValue()));
        }
        enumStrings = new Vector(obj.enumStrings.size());
        for(i = 0; i < obj.enumStrings.size(); ++i)
            enumStrings.addElement(new String((String)obj.enumStrings.elementAt(i)));

        enumStringsSizes = new Vector();
        for(i = 0; i < obj.enumStringsSizes.size(); ++i)
        {
            Integer iv = (Integer)obj.enumStringsSizes.elementAt(i);
            enumStringsSizes.addElement(new Integer(iv.intValue()));
        }
        obsoleteNames = new Vector(obj.obsoleteNames.size());
        for(i = 0; i < obj.obsoleteNames.size(); ++i)
            obsoleteNames.addElement(new String((String)obj.obsoleteNames.elementAt(i)));


        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return DBOptionsAttributes_numAdditionalAtts;
    }

    public boolean equals(DBOptionsAttributes obj)
    {
        int i;

        // Compare the elements in the types vector.
        boolean types_equal = (obj.types.size() == types.size());
        for(i = 0; (i < types.size()) && types_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer types1 = (Integer)types.elementAt(i);
            Integer types2 = (Integer)obj.types.elementAt(i);
            types_equal = types1.equals(types2);
        }
        // Compare the elements in the names vector.
        boolean names_equal = (obj.names.size() == names.size());
        for(i = 0; (i < names.size()) && names_equal; ++i)
        {
            // Make references to String from Object.
            String names1 = (String)names.elementAt(i);
            String names2 = (String)obj.names.elementAt(i);
            names_equal = names1.equals(names2);
        }
        // Compare the elements in the optBools vector.
        boolean optBools_equal = (obj.optBools.size() == optBools.size());
        for(i = 0; (i < optBools.size()) && optBools_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer optBools1 = (Integer)optBools.elementAt(i);
            Integer optBools2 = (Integer)obj.optBools.elementAt(i);
            optBools_equal = optBools1.equals(optBools2);
        }
        // Compare the elements in the optFloats vector.
        boolean optFloats_equal = (obj.optFloats.size() == optFloats.size());
        for(i = 0; (i < optFloats.size()) && optFloats_equal; ++i)
        {
            // Make references to Double from Object.
            Double optFloats1 = (Double)optFloats.elementAt(i);
            Double optFloats2 = (Double)obj.optFloats.elementAt(i);
            optFloats_equal = optFloats1.equals(optFloats2);
        }
        // Compare the elements in the optDoubles vector.
        boolean optDoubles_equal = (obj.optDoubles.size() == optDoubles.size());
        for(i = 0; (i < optDoubles.size()) && optDoubles_equal; ++i)
        {
            // Make references to Double from Object.
            Double optDoubles1 = (Double)optDoubles.elementAt(i);
            Double optDoubles2 = (Double)obj.optDoubles.elementAt(i);
            optDoubles_equal = optDoubles1.equals(optDoubles2);
        }
        // Compare the elements in the optInts vector.
        boolean optInts_equal = (obj.optInts.size() == optInts.size());
        for(i = 0; (i < optInts.size()) && optInts_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer optInts1 = (Integer)optInts.elementAt(i);
            Integer optInts2 = (Integer)obj.optInts.elementAt(i);
            optInts_equal = optInts1.equals(optInts2);
        }
        // Compare the elements in the optStrings vector.
        boolean optStrings_equal = (obj.optStrings.size() == optStrings.size());
        for(i = 0; (i < optStrings.size()) && optStrings_equal; ++i)
        {
            // Make references to String from Object.
            String optStrings1 = (String)optStrings.elementAt(i);
            String optStrings2 = (String)obj.optStrings.elementAt(i);
            optStrings_equal = optStrings1.equals(optStrings2);
        }
        // Compare the elements in the optEnums vector.
        boolean optEnums_equal = (obj.optEnums.size() == optEnums.size());
        for(i = 0; (i < optEnums.size()) && optEnums_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer optEnums1 = (Integer)optEnums.elementAt(i);
            Integer optEnums2 = (Integer)obj.optEnums.elementAt(i);
            optEnums_equal = optEnums1.equals(optEnums2);
        }
        // Compare the elements in the enumStrings vector.
        boolean enumStrings_equal = (obj.enumStrings.size() == enumStrings.size());
        for(i = 0; (i < enumStrings.size()) && enumStrings_equal; ++i)
        {
            // Make references to String from Object.
            String enumStrings1 = (String)enumStrings.elementAt(i);
            String enumStrings2 = (String)obj.enumStrings.elementAt(i);
            enumStrings_equal = enumStrings1.equals(enumStrings2);
        }
        // Compare the elements in the enumStringsSizes vector.
        boolean enumStringsSizes_equal = (obj.enumStringsSizes.size() == enumStringsSizes.size());
        for(i = 0; (i < enumStringsSizes.size()) && enumStringsSizes_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer enumStringsSizes1 = (Integer)enumStringsSizes.elementAt(i);
            Integer enumStringsSizes2 = (Integer)obj.enumStringsSizes.elementAt(i);
            enumStringsSizes_equal = enumStringsSizes1.equals(enumStringsSizes2);
        }
        // Compare the elements in the obsoleteNames vector.
        boolean obsoleteNames_equal = (obj.obsoleteNames.size() == obsoleteNames.size());
        for(i = 0; (i < obsoleteNames.size()) && obsoleteNames_equal; ++i)
        {
            // Make references to String from Object.
            String obsoleteNames1 = (String)obsoleteNames.elementAt(i);
            String obsoleteNames2 = (String)obj.obsoleteNames.elementAt(i);
            obsoleteNames_equal = obsoleteNames1.equals(obsoleteNames2);
        }
        // Create the return value
        return (types_equal &&
                names_equal &&
                optBools_equal &&
                optFloats_equal &&
                optDoubles_equal &&
                optInts_equal &&
                optStrings_equal &&
                optEnums_equal &&
                enumStrings_equal &&
                enumStringsSizes_equal &&
                obsoleteNames_equal);
    }

    // Property setting methods
    public void SetTypes(Vector types_)
    {
        types = types_;
        Select(0);
    }

    public void SetNames(Vector names_)
    {
        names = names_;
        Select(1);
    }

    public void SetOptBools(Vector optBools_)
    {
        optBools = optBools_;
        Select(2);
    }

    public void SetOptFloats(Vector optFloats_)
    {
        optFloats = optFloats_;
        Select(3);
    }

    public void SetOptDoubles(Vector optDoubles_)
    {
        optDoubles = optDoubles_;
        Select(4);
    }

    public void SetOptInts(Vector optInts_)
    {
        optInts = optInts_;
        Select(5);
    }

    public void SetOptStrings(Vector optStrings_)
    {
        optStrings = optStrings_;
        Select(6);
    }

    public void SetOptEnums(Vector optEnums_)
    {
        optEnums = optEnums_;
        Select(7);
    }

    public void SetEnumStrings(Vector enumStrings_)
    {
        enumStrings = enumStrings_;
        Select(8);
    }

    public void SetEnumStringsSizes(Vector enumStringsSizes_)
    {
        enumStringsSizes = enumStringsSizes_;
        Select(9);
    }

    public void SetObsoleteNames(Vector obsoleteNames_)
    {
        obsoleteNames = obsoleteNames_;
        Select(10);
    }

    // Property getting methods
    public Vector GetTypes() { return types; }
    public Vector GetNames() { return names; }
    public Vector GetOptBools() { return optBools; }
    public Vector GetOptFloats() { return optFloats; }
    public Vector GetOptDoubles() { return optDoubles; }
    public Vector GetOptInts() { return optInts; }
    public Vector GetOptStrings() { return optStrings; }
    public Vector GetOptEnums() { return optEnums; }
    public Vector GetEnumStrings() { return enumStrings; }
    public Vector GetEnumStringsSizes() { return enumStringsSizes; }
    public Vector GetObsoleteNames() { return obsoleteNames; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteIntVector(types);
        if(WriteSelect(1, buf))
            buf.WriteStringVector(names);
        if(WriteSelect(2, buf))
            buf.WriteIntVector(optBools);
        if(WriteSelect(3, buf))
            buf.WriteDoubleVector(optFloats);
        if(WriteSelect(4, buf))
            buf.WriteDoubleVector(optDoubles);
        if(WriteSelect(5, buf))
            buf.WriteIntVector(optInts);
        if(WriteSelect(6, buf))
            buf.WriteStringVector(optStrings);
        if(WriteSelect(7, buf))
            buf.WriteIntVector(optEnums);
        if(WriteSelect(8, buf))
            buf.WriteStringVector(enumStrings);
        if(WriteSelect(9, buf))
            buf.WriteIntVector(enumStringsSizes);
        if(WriteSelect(10, buf))
            buf.WriteStringVector(obsoleteNames);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetTypes(buf.ReadIntVector());
            break;
        case 1:
            SetNames(buf.ReadStringVector());
            break;
        case 2:
            SetOptBools(buf.ReadIntVector());
            break;
        case 3:
            SetOptFloats(buf.ReadDoubleVector());
            break;
        case 4:
            SetOptDoubles(buf.ReadDoubleVector());
            break;
        case 5:
            SetOptInts(buf.ReadIntVector());
            break;
        case 6:
            SetOptStrings(buf.ReadStringVector());
            break;
        case 7:
            SetOptEnums(buf.ReadIntVector());
            break;
        case 8:
            SetEnumStrings(buf.ReadStringVector());
            break;
        case 9:
            SetEnumStringsSizes(buf.ReadIntVector());
            break;
        case 10:
            SetObsoleteNames(buf.ReadStringVector());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + intVectorToString("types", types, indent) + "\n";
        str = str + stringVectorToString("names", names, indent) + "\n";
        str = str + intVectorToString("optBools", optBools, indent) + "\n";
        str = str + doubleVectorToString("optFloats", optFloats, indent) + "\n";
        str = str + doubleVectorToString("optDoubles", optDoubles, indent) + "\n";
        str = str + intVectorToString("optInts", optInts, indent) + "\n";
        str = str + stringVectorToString("optStrings", optStrings, indent) + "\n";
        str = str + intVectorToString("optEnums", optEnums, indent) + "\n";
        str = str + stringVectorToString("enumStrings", enumStrings, indent) + "\n";
        str = str + intVectorToString("enumStringsSizes", enumStringsSizes, indent) + "\n";
        str = str + stringVectorToString("obsoleteNames", obsoleteNames, indent) + "\n";
        return str;
    }


    // Attributes
    private Vector types; // vector of Integer objects
    private Vector names; // vector of String objects
    private Vector optBools; // vector of Integer objects
    private Vector optFloats; // vector of Double objects
    private Vector optDoubles; // vector of Double objects
    private Vector optInts; // vector of Integer objects
    private Vector optStrings; // vector of String objects
    private Vector optEnums; // vector of Integer objects
    private Vector enumStrings; // vector of String objects
    private Vector enumStringsSizes; // vector of Integer objects
    private Vector obsoleteNames; // vector of String objects
}

