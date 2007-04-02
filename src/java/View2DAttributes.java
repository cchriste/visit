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

package llnl.visit;


// ****************************************************************************
// Class: View2DAttributes
//
// Purpose:
//    This class contains the 2d view attributes.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Tue Mar 7 13:47:09 PST 2006
//
// Modifications:
//   
// ****************************************************************************

public class View2DAttributes extends AttributeSubject
{
    // Enum values
    public final static int TRISTATEMODE_ON = 0;
    public final static int TRISTATEMODE_OFF = 1;
    public final static int TRISTATEMODE_AUTO = 2;

    // Constants
    public final static double DEFAULT_FULL_FRAME_AUTO_THRESHOLD = 100.0;


    public View2DAttributes()
    {
        super(4);

        windowCoords = new double[4];
        windowCoords[0] = 0;
        windowCoords[1] = 0;
        windowCoords[2] = 1;
        windowCoords[3] = 1;
        viewportCoords = new double[4];
        viewportCoords[0] = 0.1;
        viewportCoords[1] = 0.1;
        viewportCoords[2] = 0.9;
        viewportCoords[3] = 0.9;
        fullFrameActivationMode = TRISTATEMODE_AUTO;
        fullFrameAutoThreshold = 100;
    }

    public View2DAttributes(View2DAttributes obj)
    {
        super(4);

        int i;

        windowCoords = new double[4];
        for(i = 0; i < obj.windowCoords.length; ++i)
            windowCoords[i] = obj.windowCoords[i];

        viewportCoords = new double[4];
        for(i = 0; i < obj.viewportCoords.length; ++i)
            viewportCoords[i] = obj.viewportCoords[i];

        fullFrameActivationMode = obj.fullFrameActivationMode;
        fullFrameAutoThreshold = obj.fullFrameAutoThreshold;

        SelectAll();
    }

    public boolean equals(View2DAttributes obj)
    {
        int i;

        // Compare the windowCoords arrays.
        boolean windowCoords_equal = true;
        for(i = 0; i < 4 && windowCoords_equal; ++i)
            windowCoords_equal = (windowCoords[i] == obj.windowCoords[i]);

        // Compare the viewportCoords arrays.
        boolean viewportCoords_equal = true;
        for(i = 0; i < 4 && viewportCoords_equal; ++i)
            viewportCoords_equal = (viewportCoords[i] == obj.viewportCoords[i]);

        // Create the return value
        return (windowCoords_equal &&
                viewportCoords_equal &&
                (fullFrameActivationMode == obj.fullFrameActivationMode) &&
                (fullFrameAutoThreshold == obj.fullFrameAutoThreshold));
    }

    // Property setting methods
    public void SetWindowCoords(double[] windowCoords_)
    {
        windowCoords[0] = windowCoords_[0];
        windowCoords[1] = windowCoords_[1];
        windowCoords[2] = windowCoords_[2];
        windowCoords[3] = windowCoords_[3];
        Select(0);
    }

    public void SetWindowCoords(double e0, double e1, double e2, double e3)
    {
        windowCoords[0] = e0;
        windowCoords[1] = e1;
        windowCoords[2] = e2;
        windowCoords[3] = e3;
        Select(0);
    }

    public void SetViewportCoords(double[] viewportCoords_)
    {
        viewportCoords[0] = viewportCoords_[0];
        viewportCoords[1] = viewportCoords_[1];
        viewportCoords[2] = viewportCoords_[2];
        viewportCoords[3] = viewportCoords_[3];
        Select(1);
    }

    public void SetViewportCoords(double e0, double e1, double e2, double e3)
    {
        viewportCoords[0] = e0;
        viewportCoords[1] = e1;
        viewportCoords[2] = e2;
        viewportCoords[3] = e3;
        Select(1);
    }

    public void SetFullFrameActivationMode(int fullFrameActivationMode_)
    {
        fullFrameActivationMode = fullFrameActivationMode_;
        Select(2);
    }

    public void SetFullFrameAutoThreshold(double fullFrameAutoThreshold_)
    {
        fullFrameAutoThreshold = fullFrameAutoThreshold_;
        Select(3);
    }

    // Property getting methods
    public double[] GetWindowCoords() { return windowCoords; }
    public double[] GetViewportCoords() { return viewportCoords; }
    public int      GetFullFrameActivationMode() { return fullFrameActivationMode; }
    public double   GetFullFrameAutoThreshold() { return fullFrameAutoThreshold; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteDoubleArray(windowCoords);
        if(WriteSelect(1, buf))
            buf.WriteDoubleArray(viewportCoords);
        if(WriteSelect(2, buf))
            buf.WriteInt(fullFrameActivationMode);
        if(WriteSelect(3, buf))
            buf.WriteDouble(fullFrameAutoThreshold);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetWindowCoords(buf.ReadDoubleArray());
                break;
            case 1:
                SetViewportCoords(buf.ReadDoubleArray());
                break;
            case 2:
                SetFullFrameActivationMode(buf.ReadInt());
                break;
            case 3:
                SetFullFrameAutoThreshold(buf.ReadDouble());
                break;
            }
        }
    }


    // Attributes
    private double[] windowCoords;
    private double[] viewportCoords;
    private int      fullFrameActivationMode;
    private double   fullFrameAutoThreshold;
}

