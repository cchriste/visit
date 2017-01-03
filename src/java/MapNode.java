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

package llnl.visit;

import java.util.Vector;

// ****************************************************************************
// Class: MapNode
//
// Purpose:
//   Contains a map of variant data -- a dictionary.
//
// Notes:      Maybe this would be better written using Java's HashTable...
//
// Programmer: Brad Whitlock
// Creation:   Thu Feb  2 10:26:33 PST 2012
//
// Modifications:
//   
// ****************************************************************************

public class MapNode extends Variant
{
    public MapNode()
    {
        super();
        entries = new Vector();
    }

    public MapNode(MapNode obj)
    {
        super(obj);
        entries = new Vector(obj.entries);
    }

    public MapNode(boolean value)
    {
        super();
        entries = new Vector();
        SetValue(value);
    }

    public MapNode(byte value)
    {
        super();
        entries = new Vector();
        SetValue(value);
    }

    public MapNode(int value)
    {
        super();
        entries = new Vector();
        SetValue(value);
    }

    public MapNode(long value)
    {
        super();
        entries = new Vector();
        SetValue(value);
    }

    public MapNode(float value)
    {
        super();
        entries = new Vector();
        SetValue(value);
    }

    public MapNode(double value)
    {
        super();
        entries = new Vector();
        SetValue(value);
    }

    public MapNode(String value)
    {
        super();
        entries = new Vector();
        SetValue(value);
    }

    // Write and read methods.
    public void Write(CommunicationBuffer buf)
    {
        buf.WriteInt(Type());

        if(Type() == VARIANT_EMPTY_TYPE)
        {
            buf.WriteInt(entries.size());

            for(int i = 0; i < entries.size(); ++i)
            {
                MapNodePair m = (MapNodePair)entries.elementAt(i);
                buf.WriteString(m.key);

                m.value.Write(buf);
            }
        }
        else
            super.Write(buf);
    }

    public void Read(CommunicationBuffer buf)
    {
        entries.clear();

        // Read the data type.
        dataType = buf.ReadInt();
        dataValue = null;

        if(dataType == VARIANT_EMPTY_TYPE)
        {
            int nEntries = buf.ReadInt();

            for(int i = 0; i < nEntries; ++i)
            {
                MapNodePair pair = new MapNodePair();

                pair.key = buf.ReadString();
                pair.value.Read(buf);
                
                entries.addElement(pair);
            }
        }
        else
            super.Read(buf);
    }

    public MapNode Lookup(String key)
    {
        MapNode retval = null;

        if(dataType == VARIANT_EMPTY_TYPE)
        {
            for(int i = 0; i < entries.size(); ++i)
            {
                MapNodePair pair = (MapNodePair)entries.elementAt(i);
                if(pair.key.equals(key))
                {
                    retval = pair.value;
                    break;
                }
            }
        }

        return retval;
    }

    public void SetValue(String key, MapNode value)
    {
        MapNode dest = Lookup(key);
        if(dest == null)
        {
            MapNodePair pair = new MapNodePair();
            pair.key = key; //new String(key);
            pair.value = value; //new MapNode(value);

            entries.addElement(pair);
        }
        else
            dest = value;
    }

    public String toString(String indent)
    {
        String str = new String();

        if(dataType == VARIANT_EMPTY_TYPE)
        {
            str = str + indent + "{\n";
            String indent2 = new String(indent + "    ");
            for(int i = 0; i < entries.size(); ++i)
            {
                MapNodePair pair = (MapNodePair)entries.elementAt(i);
                if(pair.value.Type() == VARIANT_EMPTY_TYPE)
                {
                    str = str + indent2 + pair.key.toString() + " = {\n";
                    str = str + pair.value.toString(indent2 + "    ") + "\n";
                    str = str + indent2 + "}\n";
                }
                else
                    str = str + indent2 + pair.key + " = " + pair.value.toString("") + "\n";
            }
            str = str + indent + "}";
        }
        else
            str = super.toString(indent);

        return str;
    }

    private Vector entries;
}
