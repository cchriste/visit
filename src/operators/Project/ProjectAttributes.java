// ***************************************************************************
//
// Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
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

// ****************************************************************************
// Class: ProjectAttributes
//
// Purpose:
//    Project data from three to two dimensions
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ProjectAttributes extends AttributeSubject implements Plugin
{
    private static int ProjectAttributes_numAdditionalAtts = 2;

    // Enum values
    public final static int PROJECTIONTYPE_ZYCARTESIAN = 0;
    public final static int PROJECTIONTYPE_XZCARTESIAN = 1;
    public final static int PROJECTIONTYPE_XYCARTESIAN = 2;
    public final static int PROJECTIONTYPE_XRCYLINDRICAL = 3;
    public final static int PROJECTIONTYPE_YRCYLINDRICAL = 4;
    public final static int PROJECTIONTYPE_ZRCYLINDRICAL = 5;

    public final static int VECTORTRANSFORMMETHOD_NONE = 0;
    public final static int VECTORTRANSFORMMETHOD_ASPOINT = 1;
    public final static int VECTORTRANSFORMMETHOD_ASDISPLACEMENT = 2;
    public final static int VECTORTRANSFORMMETHOD_ASDIRECTION = 3;


    public ProjectAttributes()
    {
        super(ProjectAttributes_numAdditionalAtts);

        projectionType = PROJECTIONTYPE_XYCARTESIAN;
        vectorTransformMethod = VECTORTRANSFORMMETHOD_ASDIRECTION;
    }

    public ProjectAttributes(int nMoreFields)
    {
        super(ProjectAttributes_numAdditionalAtts + nMoreFields);

        projectionType = PROJECTIONTYPE_XYCARTESIAN;
        vectorTransformMethod = VECTORTRANSFORMMETHOD_ASDIRECTION;
    }

    public ProjectAttributes(ProjectAttributes obj)
    {
        super(ProjectAttributes_numAdditionalAtts);

        projectionType = obj.projectionType;
        vectorTransformMethod = obj.vectorTransformMethod;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ProjectAttributes_numAdditionalAtts;
    }

    public boolean equals(ProjectAttributes obj)
    {
        // Create the return value
        return ((projectionType == obj.projectionType) &&
                (vectorTransformMethod == obj.vectorTransformMethod));
    }

    public String GetName() { return "Project"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetProjectionType(int projectionType_)
    {
        projectionType = projectionType_;
        Select(0);
    }

    public void SetVectorTransformMethod(int vectorTransformMethod_)
    {
        vectorTransformMethod = vectorTransformMethod_;
        Select(1);
    }

    // Property getting methods
    public int GetProjectionType() { return projectionType; }
    public int GetVectorTransformMethod() { return vectorTransformMethod; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(projectionType);
        if(WriteSelect(1, buf))
            buf.WriteInt(vectorTransformMethod);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetProjectionType(buf.ReadInt());
            break;
        case 1:
            SetVectorTransformMethod(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "projectionType = ";
        if(projectionType == PROJECTIONTYPE_ZYCARTESIAN)
            str = str + "PROJECTIONTYPE_ZYCARTESIAN";
        if(projectionType == PROJECTIONTYPE_XZCARTESIAN)
            str = str + "PROJECTIONTYPE_XZCARTESIAN";
        if(projectionType == PROJECTIONTYPE_XYCARTESIAN)
            str = str + "PROJECTIONTYPE_XYCARTESIAN";
        if(projectionType == PROJECTIONTYPE_XRCYLINDRICAL)
            str = str + "PROJECTIONTYPE_XRCYLINDRICAL";
        if(projectionType == PROJECTIONTYPE_YRCYLINDRICAL)
            str = str + "PROJECTIONTYPE_YRCYLINDRICAL";
        if(projectionType == PROJECTIONTYPE_ZRCYLINDRICAL)
            str = str + "PROJECTIONTYPE_ZRCYLINDRICAL";
        str = str + "\n";
        str = str + indent + "vectorTransformMethod = ";
        if(vectorTransformMethod == VECTORTRANSFORMMETHOD_NONE)
            str = str + "VECTORTRANSFORMMETHOD_NONE";
        if(vectorTransformMethod == VECTORTRANSFORMMETHOD_ASPOINT)
            str = str + "VECTORTRANSFORMMETHOD_ASPOINT";
        if(vectorTransformMethod == VECTORTRANSFORMMETHOD_ASDISPLACEMENT)
            str = str + "VECTORTRANSFORMMETHOD_ASDISPLACEMENT";
        if(vectorTransformMethod == VECTORTRANSFORMMETHOD_ASDIRECTION)
            str = str + "VECTORTRANSFORMMETHOD_ASDIRECTION";
        str = str + "\n";
        return str;
    }


    // Attributes
    private int projectionType;
    private int vectorTransformMethod;
}

