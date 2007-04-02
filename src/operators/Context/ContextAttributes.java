// ****************************************************************************
//
// Copyright (c) 2000 - 2007, The Regents of the University of California
// Produced at the Lawrence Livermore National Laboratory
// All rights reserved.
//
// This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
//    documentation and/or materials provided with the distribution.
//  - Neither the name of the UC/LLNL nor  the names of its contributors may be
//    used to  endorse or  promote products derived from  this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
// CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
// ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
// SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
// CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
// LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
// OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ****************************************************************************

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: ContextAttributes
//
// Purpose:
//    This class contains attributes for the context operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Thu Jul 31 16:09:06 PST 2003
//
// Modifications:
//   
// ****************************************************************************

public class ContextAttributes extends AttributeSubject implements Plugin
{
    // Constants
    public final static int AMOUNT_SOME = 0;
    public final static int AMOUNT_ALL = 1;


    public ContextAttributes()
    {
        super(7);

        offset = 2;
        low = 2;
        hi = 2;
        context = new String("");
        cutoff = 0;
        below = 0;
        above = 1;
    }

    public ContextAttributes(ContextAttributes obj)
    {
        super(7);

        offset = obj.offset;
        low = obj.low;
        hi = obj.hi;
        context = new String(obj.context);
        cutoff = obj.cutoff;
        below = obj.below;
        above = obj.above;

        SelectAll();
    }

    public boolean equals(ContextAttributes obj)
    {
        // Create the return value
        return ((offset == obj.offset) &&
                (low == obj.low) &&
                (hi == obj.hi) &&
                (context == obj.context) &&
                (cutoff == obj.cutoff) &&
                (below == obj.below) &&
                (above == obj.above));
    }

    public String GetName() { return "Context"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetOffset(double offset_)
    {
        offset = offset_;
        Select(0);
    }

    public void SetLow(double low_)
    {
        low = low_;
        Select(1);
    }

    public void SetHi(double hi_)
    {
        hi = hi_;
        Select(2);
    }

    public void SetContext(String context_)
    {
        context = context_;
        Select(3);
    }

    public void SetCutoff(double cutoff_)
    {
        cutoff = cutoff_;
        Select(4);
    }

    public void SetBelow(double below_)
    {
        below = below_;
        Select(5);
    }

    public void SetAbove(double above_)
    {
        above = above_;
        Select(6);
    }

    // Property getting methods
    public double GetOffset() { return offset; }
    public double GetLow() { return low; }
    public double GetHi() { return hi; }
    public String GetContext() { return context; }
    public double GetCutoff() { return cutoff; }
    public double GetBelow() { return below; }
    public double GetAbove() { return above; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDouble(offset);
        if(WriteSelect(1, buf))
            buf.WriteDouble(low);
        if(WriteSelect(2, buf))
            buf.WriteDouble(hi);
        if(WriteSelect(3, buf))
            buf.WriteString(context);
        if(WriteSelect(4, buf))
            buf.WriteDouble(cutoff);
        if(WriteSelect(5, buf))
            buf.WriteDouble(below);
        if(WriteSelect(6, buf))
            buf.WriteDouble(above);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetOffset(buf.ReadDouble());
                break;
            case 1:
                SetLow(buf.ReadDouble());
                break;
            case 2:
                SetHi(buf.ReadDouble());
                break;
            case 3:
                SetContext(buf.ReadString());
                break;
            case 4:
                SetCutoff(buf.ReadDouble());
                break;
            case 5:
                SetBelow(buf.ReadDouble());
                break;
            case 6:
                SetAbove(buf.ReadDouble());
                break;
            }
        }
    }


    // Attributes
    private double offset;
    private double low;
    private double hi;
    private String context;
    private double cutoff;
    private double below;
    private double above;
}

