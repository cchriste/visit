// ***************************************************************************
//
// Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
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

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;
import java.util.Vector;
import java.lang.Integer;
import java.lang.Double;
import java.lang.Byte;

// ****************************************************************************
// Class: ModelFitAtts
//
// Purpose:
//    This file contains attributes for the ModelFit operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ModelFitAtts extends AttributeSubject implements Plugin
{
    private static int ModelFitAtts_numAdditionalAtts = 11;

    // Enum values
    public final static int SPACES_VARIABLE = 0;
    public final static int SPACES_NORMAL = 1;
    public final static int SPACES_LOG = 2;
    public final static int SPACES_PROBABILITY = 3;

    public final static int STATS_AVERAGE = 0;
    public final static int STATS_MIN = 1;
    public final static int STATS_MAX = 2;
    public final static int STATS_NONE = 3;

    public final static int DISTANCES_EUCLIDEAN = 0;
    public final static int DISTANCES_MANHATTAN = 1;
    public final static int DISTANCES_MAXIMUM = 2;


    public ModelFitAtts()
    {
        super(ModelFitAtts_numAdditionalAtts);

        Vars = new Vector();
        numVars = new Vector();
        Tuples = new Vector();
        StatTuples = new Vector();
        numTups = new Vector();
        thold = new Vector();
        selectionType = new Vector();
        distanceType = new Vector();
        inputSpace = new Vector();
        modelNames = new Vector();
        modelNums = new Vector();
    }

    public ModelFitAtts(int nMoreFields)
    {
        super(ModelFitAtts_numAdditionalAtts + nMoreFields);

        Vars = new Vector();
        numVars = new Vector();
        Tuples = new Vector();
        StatTuples = new Vector();
        numTups = new Vector();
        thold = new Vector();
        selectionType = new Vector();
        distanceType = new Vector();
        inputSpace = new Vector();
        modelNames = new Vector();
        modelNums = new Vector();
    }

    public ModelFitAtts(ModelFitAtts obj)
    {
        super(obj);

        int i;

        Vars = new Vector(obj.Vars.size());
        for(i = 0; i < obj.Vars.size(); ++i)
            Vars.addElement(new String((String)obj.Vars.elementAt(i)));

        numVars = new Vector();
        for(i = 0; i < obj.numVars.size(); ++i)
        {
            Integer iv = (Integer)obj.numVars.elementAt(i);
            numVars.addElement(new Integer(iv.intValue()));
        }
        Tuples = new Vector(obj.Tuples.size());
        for(i = 0; i < obj.Tuples.size(); ++i)
        {
            Double dv = (Double)obj.Tuples.elementAt(i);
            Tuples.addElement(new Double(dv.doubleValue()));
        }

        StatTuples = new Vector(obj.StatTuples.size());
        for(i = 0; i < obj.StatTuples.size(); ++i)
        {
            Byte bv = (Byte)obj.StatTuples.elementAt(i);
            StatTuples.addElement(new Byte(bv.byteValue()));
        }

        numTups = new Vector();
        for(i = 0; i < obj.numTups.size(); ++i)
        {
            Integer iv = (Integer)obj.numTups.elementAt(i);
            numTups.addElement(new Integer(iv.intValue()));
        }
        thold = new Vector(obj.thold.size());
        for(i = 0; i < obj.thold.size(); ++i)
        {
            Double dv = (Double)obj.thold.elementAt(i);
            thold.addElement(new Double(dv.doubleValue()));
        }

        selectionType = new Vector();
        for(i = 0; i < obj.selectionType.size(); ++i)
        {
            Integer iv = (Integer)obj.selectionType.elementAt(i);
            selectionType.addElement(new Integer(iv.intValue()));
        }
        distanceType = new Vector();
        for(i = 0; i < obj.distanceType.size(); ++i)
        {
            Integer iv = (Integer)obj.distanceType.elementAt(i);
            distanceType.addElement(new Integer(iv.intValue()));
        }
        inputSpace = new Vector();
        for(i = 0; i < obj.inputSpace.size(); ++i)
        {
            Integer iv = (Integer)obj.inputSpace.elementAt(i);
            inputSpace.addElement(new Integer(iv.intValue()));
        }
        modelNames = new Vector(obj.modelNames.size());
        for(i = 0; i < obj.modelNames.size(); ++i)
            modelNames.addElement(new String((String)obj.modelNames.elementAt(i)));

        modelNums = new Vector();
        for(i = 0; i < obj.modelNums.size(); ++i)
        {
            Integer iv = (Integer)obj.modelNums.elementAt(i);
            modelNums.addElement(new Integer(iv.intValue()));
        }

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ModelFitAtts_numAdditionalAtts;
    }

    public boolean equals(ModelFitAtts obj)
    {
        int i;

        // Compare the elements in the Vars vector.
        boolean Vars_equal = (obj.Vars.size() == Vars.size());
        for(i = 0; (i < Vars.size()) && Vars_equal; ++i)
        {
            // Make references to String from Object.
            String Vars1 = (String)Vars.elementAt(i);
            String Vars2 = (String)obj.Vars.elementAt(i);
            Vars_equal = Vars1.equals(Vars2);
        }
        // Compare the elements in the numVars vector.
        boolean numVars_equal = (obj.numVars.size() == numVars.size());
        for(i = 0; (i < numVars.size()) && numVars_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer numVars1 = (Integer)numVars.elementAt(i);
            Integer numVars2 = (Integer)obj.numVars.elementAt(i);
            numVars_equal = numVars1.equals(numVars2);
        }
        // Compare the elements in the Tuples vector.
        boolean Tuples_equal = (obj.Tuples.size() == Tuples.size());
        for(i = 0; (i < Tuples.size()) && Tuples_equal; ++i)
        {
            // Make references to Double from Object.
            Double Tuples1 = (Double)Tuples.elementAt(i);
            Double Tuples2 = (Double)obj.Tuples.elementAt(i);
            Tuples_equal = Tuples1.equals(Tuples2);
        }
        // Compare the elements in the StatTuples vector.
        boolean StatTuples_equal = (obj.StatTuples.size() == StatTuples.size());
        for(i = 0; (i < StatTuples.size()) && StatTuples_equal; ++i)
        {
            // Make references to Byte from Object.
            Byte StatTuples1 = (Byte)StatTuples.elementAt(i);
            Byte StatTuples2 = (Byte)obj.StatTuples.elementAt(i);
            StatTuples_equal = StatTuples1.equals(StatTuples2);
        }
        // Compare the elements in the numTups vector.
        boolean numTups_equal = (obj.numTups.size() == numTups.size());
        for(i = 0; (i < numTups.size()) && numTups_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer numTups1 = (Integer)numTups.elementAt(i);
            Integer numTups2 = (Integer)obj.numTups.elementAt(i);
            numTups_equal = numTups1.equals(numTups2);
        }
        // Compare the elements in the thold vector.
        boolean thold_equal = (obj.thold.size() == thold.size());
        for(i = 0; (i < thold.size()) && thold_equal; ++i)
        {
            // Make references to Double from Object.
            Double thold1 = (Double)thold.elementAt(i);
            Double thold2 = (Double)obj.thold.elementAt(i);
            thold_equal = thold1.equals(thold2);
        }
        // Compare the elements in the selectionType vector.
        boolean selectionType_equal = (obj.selectionType.size() == selectionType.size());
        for(i = 0; (i < selectionType.size()) && selectionType_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer selectionType1 = (Integer)selectionType.elementAt(i);
            Integer selectionType2 = (Integer)obj.selectionType.elementAt(i);
            selectionType_equal = selectionType1.equals(selectionType2);
        }
        // Compare the elements in the distanceType vector.
        boolean distanceType_equal = (obj.distanceType.size() == distanceType.size());
        for(i = 0; (i < distanceType.size()) && distanceType_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer distanceType1 = (Integer)distanceType.elementAt(i);
            Integer distanceType2 = (Integer)obj.distanceType.elementAt(i);
            distanceType_equal = distanceType1.equals(distanceType2);
        }
        // Compare the elements in the inputSpace vector.
        boolean inputSpace_equal = (obj.inputSpace.size() == inputSpace.size());
        for(i = 0; (i < inputSpace.size()) && inputSpace_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer inputSpace1 = (Integer)inputSpace.elementAt(i);
            Integer inputSpace2 = (Integer)obj.inputSpace.elementAt(i);
            inputSpace_equal = inputSpace1.equals(inputSpace2);
        }
        // Compare the elements in the modelNames vector.
        boolean modelNames_equal = (obj.modelNames.size() == modelNames.size());
        for(i = 0; (i < modelNames.size()) && modelNames_equal; ++i)
        {
            // Make references to String from Object.
            String modelNames1 = (String)modelNames.elementAt(i);
            String modelNames2 = (String)obj.modelNames.elementAt(i);
            modelNames_equal = modelNames1.equals(modelNames2);
        }
        // Compare the elements in the modelNums vector.
        boolean modelNums_equal = (obj.modelNums.size() == modelNums.size());
        for(i = 0; (i < modelNums.size()) && modelNums_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer modelNums1 = (Integer)modelNums.elementAt(i);
            Integer modelNums2 = (Integer)obj.modelNums.elementAt(i);
            modelNums_equal = modelNums1.equals(modelNums2);
        }
        // Create the return value
        return (Vars_equal &&
                numVars_equal &&
                Tuples_equal &&
                StatTuples_equal &&
                numTups_equal &&
                thold_equal &&
                selectionType_equal &&
                distanceType_equal &&
                inputSpace_equal &&
                modelNames_equal &&
                modelNums_equal);
    }

    public String GetName() { return "ModelFit"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetVars(Vector Vars_)
    {
        Vars = Vars_;
        Select(0);
    }

    public void SetNumVars(Vector numVars_)
    {
        numVars = numVars_;
        Select(1);
    }

    public void SetTuples(Vector Tuples_)
    {
        Tuples = Tuples_;
        Select(2);
    }

    public void SetStatTuples(Vector StatTuples_)
    {
        StatTuples = StatTuples_;
        Select(3);
    }

    public void SetNumTups(Vector numTups_)
    {
        numTups = numTups_;
        Select(4);
    }

    public void SetThold(Vector thold_)
    {
        thold = thold_;
        Select(5);
    }

    public void SetSelectionType(Vector selectionType_)
    {
        selectionType = selectionType_;
        Select(6);
    }

    public void SetDistanceType(Vector distanceType_)
    {
        distanceType = distanceType_;
        Select(7);
    }

    public void SetInputSpace(Vector inputSpace_)
    {
        inputSpace = inputSpace_;
        Select(8);
    }

    public void SetModelNames(Vector modelNames_)
    {
        modelNames = modelNames_;
        Select(9);
    }

    public void SetModelNums(Vector modelNums_)
    {
        modelNums = modelNums_;
        Select(10);
    }

    // Property getting methods
    public Vector GetVars() { return Vars; }
    public Vector GetNumVars() { return numVars; }
    public Vector GetTuples() { return Tuples; }
    public Vector GetStatTuples() { return StatTuples; }
    public Vector GetNumTups() { return numTups; }
    public Vector GetThold() { return thold; }
    public Vector GetSelectionType() { return selectionType; }
    public Vector GetDistanceType() { return distanceType; }
    public Vector GetInputSpace() { return inputSpace; }
    public Vector GetModelNames() { return modelNames; }
    public Vector GetModelNums() { return modelNums; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteStringVector(Vars);
        if(WriteSelect(1, buf))
            buf.WriteIntVector(numVars);
        if(WriteSelect(2, buf))
            buf.WriteDoubleVector(Tuples);
        if(WriteSelect(3, buf))
            buf.WriteByteVector(StatTuples);
        if(WriteSelect(4, buf))
            buf.WriteIntVector(numTups);
        if(WriteSelect(5, buf))
            buf.WriteDoubleVector(thold);
        if(WriteSelect(6, buf))
            buf.WriteIntVector(selectionType);
        if(WriteSelect(7, buf))
            buf.WriteIntVector(distanceType);
        if(WriteSelect(8, buf))
            buf.WriteIntVector(inputSpace);
        if(WriteSelect(9, buf))
            buf.WriteStringVector(modelNames);
        if(WriteSelect(10, buf))
            buf.WriteIntVector(modelNums);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetVars(buf.ReadStringVector());
            break;
        case 1:
            SetNumVars(buf.ReadIntVector());
            break;
        case 2:
            SetTuples(buf.ReadDoubleVector());
            break;
        case 3:
            SetStatTuples(buf.ReadByteVector());
            break;
        case 4:
            SetNumTups(buf.ReadIntVector());
            break;
        case 5:
            SetThold(buf.ReadDoubleVector());
            break;
        case 6:
            SetSelectionType(buf.ReadIntVector());
            break;
        case 7:
            SetDistanceType(buf.ReadIntVector());
            break;
        case 8:
            SetInputSpace(buf.ReadIntVector());
            break;
        case 9:
            SetModelNames(buf.ReadStringVector());
            break;
        case 10:
            SetModelNums(buf.ReadIntVector());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringVectorToString("Vars", Vars, indent) + "\n";
        str = str + intVectorToString("numVars", numVars, indent) + "\n";
        str = str + doubleVectorToString("Tuples", Tuples, indent) + "\n";
        str = str + ucharVectorToString("StatTuples", StatTuples, indent) + "\n";
        str = str + intVectorToString("numTups", numTups, indent) + "\n";
        str = str + doubleVectorToString("thold", thold, indent) + "\n";
        str = str + intVectorToString("selectionType", selectionType, indent) + "\n";
        str = str + intVectorToString("distanceType", distanceType, indent) + "\n";
        str = str + intVectorToString("inputSpace", inputSpace, indent) + "\n";
        str = str + stringVectorToString("modelNames", modelNames, indent) + "\n";
        str = str + intVectorToString("modelNums", modelNums, indent) + "\n";
        return str;
    }


    // Attributes
    private Vector Vars; // vector of String objects
    private Vector numVars; // vector of Integer objects
    private Vector Tuples; // vector of Double objects
    private Vector StatTuples; // vector of Byte objects
    private Vector numTups; // vector of Integer objects
    private Vector thold; // vector of Double objects
    private Vector selectionType; // vector of Integer objects
    private Vector distanceType; // vector of Integer objects
    private Vector inputSpace; // vector of Integer objects
    private Vector modelNames; // vector of String objects
    private Vector modelNums; // vector of Integer objects
}

