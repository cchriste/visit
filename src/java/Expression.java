// ***************************************************************************
//
// Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400142
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


// ****************************************************************************
// Class: Expression
//
// Purpose:
//    This class contains an expression.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 14 17:54:52 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class Expression extends AttributeSubject
{
    // Enum values
    public final static int EXPRTYPE_UNKNOWN = 0;
    public final static int EXPRTYPE_SCALARMESHVAR = 1;
    public final static int EXPRTYPE_VECTORMESHVAR = 2;
    public final static int EXPRTYPE_TENSORMESHVAR = 3;
    public final static int EXPRTYPE_SYMMETRICTENSORMESHVAR = 4;
    public final static int EXPRTYPE_ARRAYMESHVAR = 5;
    public final static int EXPRTYPE_CURVEMESHVAR = 6;
    public final static int EXPRTYPE_MESH = 7;
    public final static int EXPRTYPE_MATERIAL = 8;
    public final static int EXPRTYPE_SPECIES = 9;


    public Expression()
    {
        super(7);

        name = new String("notset");
        definition = new String("notset");
        hidden = false;
        type = EXPRTYPE_SCALARMESHVAR;
        fromDB = false;
        dbName = new String("__none__");
        autoExpression = false;
    }

    public Expression(Expression obj)
    {
        super(7);

        name = new String(obj.name);
        definition = new String(obj.definition);
        hidden = obj.hidden;
        type = obj.type;
        fromDB = obj.fromDB;
        dbName = new String(obj.dbName);
        autoExpression = obj.autoExpression;

        SelectAll();
    }

    public boolean equals(Expression obj)
    {
        // Create the return value
        return ((name == obj.name) &&
                (definition == obj.definition) &&
                (hidden == obj.hidden) &&
                (type == obj.type) &&
                (fromDB == obj.fromDB) &&
                (dbName == obj.dbName) &&
                (autoExpression == obj.autoExpression));
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetDefinition(String definition_)
    {
        definition = definition_;
        Select(1);
    }

    public void SetHidden(boolean hidden_)
    {
        hidden = hidden_;
        Select(2);
    }

    public void SetType(int type_)
    {
        type = type_;
        Select(3);
    }

    public void SetFromDB(boolean fromDB_)
    {
        fromDB = fromDB_;
        Select(4);
    }

    public void SetDbName(String dbName_)
    {
        dbName = dbName_;
        Select(5);
    }

    public void SetAutoExpression(boolean autoExpression_)
    {
        autoExpression = autoExpression_;
        Select(6);
    }

    // Property getting methods
    public String  GetName() { return name; }
    public String  GetDefinition() { return definition; }
    public boolean GetHidden() { return hidden; }
    public int     GetType() { return type; }
    public boolean GetFromDB() { return fromDB; }
    public String  GetDbName() { return dbName; }
    public boolean GetAutoExpression() { return autoExpression; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteString(definition);
        if(WriteSelect(2, buf))
            buf.WriteBool(hidden);
        if(WriteSelect(3, buf))
            buf.WriteInt(type);
        if(WriteSelect(4, buf))
            buf.WriteBool(fromDB);
        if(WriteSelect(5, buf))
            buf.WriteString(dbName);
        if(WriteSelect(6, buf))
            buf.WriteBool(autoExpression);
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
                SetDefinition(buf.ReadString());
                break;
            case 2:
                SetHidden(buf.ReadBool());
                break;
            case 3:
                SetType(buf.ReadInt());
                break;
            case 4:
                SetFromDB(buf.ReadBool());
                break;
            case 5:
                SetDbName(buf.ReadString());
                break;
            case 6:
                SetAutoExpression(buf.ReadBool());
                break;
            }
        }
    }


    // Attributes
    private String  name;
    private String  definition;
    private boolean hidden;
    private int     type;
    private boolean fromDB;
    private String  dbName;
    private boolean autoExpression;
}

