// ***************************************************************************
//
// Copyright (c) 2000 - 2012, Lawrence Livermore National Security, LLC
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
import java.lang.Integer;

// ****************************************************************************
// Class: WindowInformation
//
// Purpose:
//    This class contains the attributes that tell the state of a viewer window.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class WindowInformation extends AttributeSubject
{
    private static int WindowInformation_numAdditionalAtts = 27;

    public WindowInformation()
    {
        super(WindowInformation_numAdditionalAtts);

        activeSource = new String("");
        activeTimeSlider = -1;
        timeSliders = new Vector();
        timeSliderCurrentStates = new Vector();
        animationMode = 2;
        interactionMode = 0;
        toolUpdateMode = 1;
        boundingBoxNavigate = true;
        spin = false;
        fullFrame = false;
        perspective = true;
        maintainView = false;
        lockView = false;
        lockTools = false;
        lockTime = false;
        viewExtentsType = 0;
        viewDimension = 2;
        viewKeyframes = new Vector();
        cameraViewMode = false;
        usingScalableRendering = false;
        lastRenderMin = 0f;
        lastRenderAvg = 0f;
        lastRenderMax = 0f;
        numPrimitives = 0;
        extents = new double[6];
        extents[0] = 0;
        extents[1] = 0;
        extents[2] = 0;
        extents[3] = 0;
        extents[4] = 0;
        extents[5] = 0;
        windowSize = new int[2];
        windowSize[0] = 0;
        windowSize[1] = 0;
        winMode = 0;
    }

    public WindowInformation(int nMoreFields)
    {
        super(WindowInformation_numAdditionalAtts + nMoreFields);

        activeSource = new String("");
        activeTimeSlider = -1;
        timeSliders = new Vector();
        timeSliderCurrentStates = new Vector();
        animationMode = 2;
        interactionMode = 0;
        toolUpdateMode = 1;
        boundingBoxNavigate = true;
        spin = false;
        fullFrame = false;
        perspective = true;
        maintainView = false;
        lockView = false;
        lockTools = false;
        lockTime = false;
        viewExtentsType = 0;
        viewDimension = 2;
        viewKeyframes = new Vector();
        cameraViewMode = false;
        usingScalableRendering = false;
        lastRenderMin = 0f;
        lastRenderAvg = 0f;
        lastRenderMax = 0f;
        numPrimitives = 0;
        extents = new double[6];
        extents[0] = 0;
        extents[1] = 0;
        extents[2] = 0;
        extents[3] = 0;
        extents[4] = 0;
        extents[5] = 0;
        windowSize = new int[2];
        windowSize[0] = 0;
        windowSize[1] = 0;
        winMode = 0;
    }

    public WindowInformation(WindowInformation obj)
    {
        super(WindowInformation_numAdditionalAtts);

        int i;

        activeSource = new String(obj.activeSource);
        activeTimeSlider = obj.activeTimeSlider;
        timeSliders = new Vector(obj.timeSliders.size());
        for(i = 0; i < obj.timeSliders.size(); ++i)
            timeSliders.addElement(new String((String)obj.timeSliders.elementAt(i)));

        timeSliderCurrentStates = new Vector();
        for(i = 0; i < obj.timeSliderCurrentStates.size(); ++i)
        {
            Integer iv = (Integer)obj.timeSliderCurrentStates.elementAt(i);
            timeSliderCurrentStates.addElement(new Integer(iv.intValue()));
        }
        animationMode = obj.animationMode;
        interactionMode = obj.interactionMode;
        toolUpdateMode = obj.toolUpdateMode;
        boundingBoxNavigate = obj.boundingBoxNavigate;
        spin = obj.spin;
        fullFrame = obj.fullFrame;
        perspective = obj.perspective;
        maintainView = obj.maintainView;
        lockView = obj.lockView;
        lockTools = obj.lockTools;
        lockTime = obj.lockTime;
        viewExtentsType = obj.viewExtentsType;
        viewDimension = obj.viewDimension;
        viewKeyframes = new Vector();
        for(i = 0; i < obj.viewKeyframes.size(); ++i)
        {
            Integer iv = (Integer)obj.viewKeyframes.elementAt(i);
            viewKeyframes.addElement(new Integer(iv.intValue()));
        }
        cameraViewMode = obj.cameraViewMode;
        usingScalableRendering = obj.usingScalableRendering;
        lastRenderMin = obj.lastRenderMin;
        lastRenderAvg = obj.lastRenderAvg;
        lastRenderMax = obj.lastRenderMax;
        numPrimitives = obj.numPrimitives;
        extents = new double[6];
        for(i = 0; i < obj.extents.length; ++i)
            extents[i] = obj.extents[i];

        windowSize = new int[2];
        windowSize[0] = obj.windowSize[0];
        windowSize[1] = obj.windowSize[1];

        winMode = obj.winMode;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return WindowInformation_numAdditionalAtts;
    }

    public boolean equals(WindowInformation obj)
    {
        int i;

        // Compare the elements in the timeSliders vector.
        boolean timeSliders_equal = (obj.timeSliders.size() == timeSliders.size());
        for(i = 0; (i < timeSliders.size()) && timeSliders_equal; ++i)
        {
            // Make references to String from Object.
            String timeSliders1 = (String)timeSliders.elementAt(i);
            String timeSliders2 = (String)obj.timeSliders.elementAt(i);
            timeSliders_equal = timeSliders1.equals(timeSliders2);
        }
        // Compare the elements in the timeSliderCurrentStates vector.
        boolean timeSliderCurrentStates_equal = (obj.timeSliderCurrentStates.size() == timeSliderCurrentStates.size());
        for(i = 0; (i < timeSliderCurrentStates.size()) && timeSliderCurrentStates_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer timeSliderCurrentStates1 = (Integer)timeSliderCurrentStates.elementAt(i);
            Integer timeSliderCurrentStates2 = (Integer)obj.timeSliderCurrentStates.elementAt(i);
            timeSliderCurrentStates_equal = timeSliderCurrentStates1.equals(timeSliderCurrentStates2);
        }
        // Compare the elements in the viewKeyframes vector.
        boolean viewKeyframes_equal = (obj.viewKeyframes.size() == viewKeyframes.size());
        for(i = 0; (i < viewKeyframes.size()) && viewKeyframes_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer viewKeyframes1 = (Integer)viewKeyframes.elementAt(i);
            Integer viewKeyframes2 = (Integer)obj.viewKeyframes.elementAt(i);
            viewKeyframes_equal = viewKeyframes1.equals(viewKeyframes2);
        }
        // Compare the extents arrays.
        boolean extents_equal = true;
        for(i = 0; i < 6 && extents_equal; ++i)
            extents_equal = (extents[i] == obj.extents[i]);

        // Compare the windowSize arrays.
        boolean windowSize_equal = true;
        for(i = 0; i < 2 && windowSize_equal; ++i)
            windowSize_equal = (windowSize[i] == obj.windowSize[i]);

        // Create the return value
        return ((activeSource.equals(obj.activeSource)) &&
                (activeTimeSlider == obj.activeTimeSlider) &&
                timeSliders_equal &&
                timeSliderCurrentStates_equal &&
                (animationMode == obj.animationMode) &&
                (interactionMode == obj.interactionMode) &&
                (toolUpdateMode == obj.toolUpdateMode) &&
                (boundingBoxNavigate == obj.boundingBoxNavigate) &&
                (spin == obj.spin) &&
                (fullFrame == obj.fullFrame) &&
                (perspective == obj.perspective) &&
                (maintainView == obj.maintainView) &&
                (lockView == obj.lockView) &&
                (lockTools == obj.lockTools) &&
                (lockTime == obj.lockTime) &&
                (viewExtentsType == obj.viewExtentsType) &&
                (viewDimension == obj.viewDimension) &&
                viewKeyframes_equal &&
                (cameraViewMode == obj.cameraViewMode) &&
                (usingScalableRendering == obj.usingScalableRendering) &&
                (lastRenderMin == obj.lastRenderMin) &&
                (lastRenderAvg == obj.lastRenderAvg) &&
                (lastRenderMax == obj.lastRenderMax) &&
                (numPrimitives == obj.numPrimitives) &&
                extents_equal &&
                windowSize_equal &&
                (winMode == obj.winMode));
    }

    // Property setting methods
    public void SetActiveSource(String activeSource_)
    {
        activeSource = activeSource_;
        Select(0);
    }

    public void SetActiveTimeSlider(int activeTimeSlider_)
    {
        activeTimeSlider = activeTimeSlider_;
        Select(1);
    }

    public void SetTimeSliders(Vector timeSliders_)
    {
        timeSliders = timeSliders_;
        Select(2);
    }

    public void SetTimeSliderCurrentStates(Vector timeSliderCurrentStates_)
    {
        timeSliderCurrentStates = timeSliderCurrentStates_;
        Select(3);
    }

    public void SetAnimationMode(int animationMode_)
    {
        animationMode = animationMode_;
        Select(4);
    }

    public void SetInteractionMode(int interactionMode_)
    {
        interactionMode = interactionMode_;
        Select(5);
    }

    public void SetToolUpdateMode(int toolUpdateMode_)
    {
        toolUpdateMode = toolUpdateMode_;
        Select(6);
    }

    public void SetBoundingBoxNavigate(boolean boundingBoxNavigate_)
    {
        boundingBoxNavigate = boundingBoxNavigate_;
        Select(7);
    }

    public void SetSpin(boolean spin_)
    {
        spin = spin_;
        Select(8);
    }

    public void SetFullFrame(boolean fullFrame_)
    {
        fullFrame = fullFrame_;
        Select(9);
    }

    public void SetPerspective(boolean perspective_)
    {
        perspective = perspective_;
        Select(10);
    }

    public void SetMaintainView(boolean maintainView_)
    {
        maintainView = maintainView_;
        Select(11);
    }

    public void SetLockView(boolean lockView_)
    {
        lockView = lockView_;
        Select(12);
    }

    public void SetLockTools(boolean lockTools_)
    {
        lockTools = lockTools_;
        Select(13);
    }

    public void SetLockTime(boolean lockTime_)
    {
        lockTime = lockTime_;
        Select(14);
    }

    public void SetViewExtentsType(int viewExtentsType_)
    {
        viewExtentsType = viewExtentsType_;
        Select(15);
    }

    public void SetViewDimension(int viewDimension_)
    {
        viewDimension = viewDimension_;
        Select(16);
    }

    public void SetViewKeyframes(Vector viewKeyframes_)
    {
        viewKeyframes = viewKeyframes_;
        Select(17);
    }

    public void SetCameraViewMode(boolean cameraViewMode_)
    {
        cameraViewMode = cameraViewMode_;
        Select(18);
    }

    public void SetUsingScalableRendering(boolean usingScalableRendering_)
    {
        usingScalableRendering = usingScalableRendering_;
        Select(19);
    }

    public void SetLastRenderMin(float lastRenderMin_)
    {
        lastRenderMin = lastRenderMin_;
        Select(20);
    }

    public void SetLastRenderAvg(float lastRenderAvg_)
    {
        lastRenderAvg = lastRenderAvg_;
        Select(21);
    }

    public void SetLastRenderMax(float lastRenderMax_)
    {
        lastRenderMax = lastRenderMax_;
        Select(22);
    }

    public void SetNumPrimitives(int numPrimitives_)
    {
        numPrimitives = numPrimitives_;
        Select(23);
    }

    public void SetExtents(double[] extents_)
    {
        for(int i = 0; i < 6; ++i)
             extents[i] = extents_[i];
        Select(24);
    }

    public void SetWindowSize(int[] windowSize_)
    {
        windowSize[0] = windowSize_[0];
        windowSize[1] = windowSize_[1];
        Select(25);
    }

    public void SetWindowSize(int e0, int e1)
    {
        windowSize[0] = e0;
        windowSize[1] = e1;
        Select(25);
    }

    public void SetWinMode(int winMode_)
    {
        winMode = winMode_;
        Select(26);
    }

    // Property getting methods
    public String   GetActiveSource() { return activeSource; }
    public int      GetActiveTimeSlider() { return activeTimeSlider; }
    public Vector   GetTimeSliders() { return timeSliders; }
    public Vector   GetTimeSliderCurrentStates() { return timeSliderCurrentStates; }
    public int      GetAnimationMode() { return animationMode; }
    public int      GetInteractionMode() { return interactionMode; }
    public int      GetToolUpdateMode() { return toolUpdateMode; }
    public boolean  GetBoundingBoxNavigate() { return boundingBoxNavigate; }
    public boolean  GetSpin() { return spin; }
    public boolean  GetFullFrame() { return fullFrame; }
    public boolean  GetPerspective() { return perspective; }
    public boolean  GetMaintainView() { return maintainView; }
    public boolean  GetLockView() { return lockView; }
    public boolean  GetLockTools() { return lockTools; }
    public boolean  GetLockTime() { return lockTime; }
    public int      GetViewExtentsType() { return viewExtentsType; }
    public int      GetViewDimension() { return viewDimension; }
    public Vector   GetViewKeyframes() { return viewKeyframes; }
    public boolean  GetCameraViewMode() { return cameraViewMode; }
    public boolean  GetUsingScalableRendering() { return usingScalableRendering; }
    public float    GetLastRenderMin() { return lastRenderMin; }
    public float    GetLastRenderAvg() { return lastRenderAvg; }
    public float    GetLastRenderMax() { return lastRenderMax; }
    public int      GetNumPrimitives() { return numPrimitives; }
    public double[] GetExtents() { return extents; }
    public int[]    GetWindowSize() { return windowSize; }
    public int      GetWinMode() { return winMode; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(activeSource);
        if(WriteSelect(1, buf))
            buf.WriteInt(activeTimeSlider);
        if(WriteSelect(2, buf))
            buf.WriteStringVector(timeSliders);
        if(WriteSelect(3, buf))
            buf.WriteIntVector(timeSliderCurrentStates);
        if(WriteSelect(4, buf))
            buf.WriteInt(animationMode);
        if(WriteSelect(5, buf))
            buf.WriteInt(interactionMode);
        if(WriteSelect(6, buf))
            buf.WriteInt(toolUpdateMode);
        if(WriteSelect(7, buf))
            buf.WriteBool(boundingBoxNavigate);
        if(WriteSelect(8, buf))
            buf.WriteBool(spin);
        if(WriteSelect(9, buf))
            buf.WriteBool(fullFrame);
        if(WriteSelect(10, buf))
            buf.WriteBool(perspective);
        if(WriteSelect(11, buf))
            buf.WriteBool(maintainView);
        if(WriteSelect(12, buf))
            buf.WriteBool(lockView);
        if(WriteSelect(13, buf))
            buf.WriteBool(lockTools);
        if(WriteSelect(14, buf))
            buf.WriteBool(lockTime);
        if(WriteSelect(15, buf))
            buf.WriteInt(viewExtentsType);
        if(WriteSelect(16, buf))
            buf.WriteInt(viewDimension);
        if(WriteSelect(17, buf))
            buf.WriteIntVector(viewKeyframes);
        if(WriteSelect(18, buf))
            buf.WriteBool(cameraViewMode);
        if(WriteSelect(19, buf))
            buf.WriteBool(usingScalableRendering);
        if(WriteSelect(20, buf))
            buf.WriteFloat(lastRenderMin);
        if(WriteSelect(21, buf))
            buf.WriteFloat(lastRenderAvg);
        if(WriteSelect(22, buf))
            buf.WriteFloat(lastRenderMax);
        if(WriteSelect(23, buf))
            buf.WriteInt(numPrimitives);
        if(WriteSelect(24, buf))
            buf.WriteDoubleArray(extents);
        if(WriteSelect(25, buf))
            buf.WriteIntArray(windowSize);
        if(WriteSelect(26, buf))
            buf.WriteInt(winMode);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetActiveSource(buf.ReadString());
            break;
        case 1:
            SetActiveTimeSlider(buf.ReadInt());
            break;
        case 2:
            SetTimeSliders(buf.ReadStringVector());
            break;
        case 3:
            SetTimeSliderCurrentStates(buf.ReadIntVector());
            break;
        case 4:
            SetAnimationMode(buf.ReadInt());
            break;
        case 5:
            SetInteractionMode(buf.ReadInt());
            break;
        case 6:
            SetToolUpdateMode(buf.ReadInt());
            break;
        case 7:
            SetBoundingBoxNavigate(buf.ReadBool());
            break;
        case 8:
            SetSpin(buf.ReadBool());
            break;
        case 9:
            SetFullFrame(buf.ReadBool());
            break;
        case 10:
            SetPerspective(buf.ReadBool());
            break;
        case 11:
            SetMaintainView(buf.ReadBool());
            break;
        case 12:
            SetLockView(buf.ReadBool());
            break;
        case 13:
            SetLockTools(buf.ReadBool());
            break;
        case 14:
            SetLockTime(buf.ReadBool());
            break;
        case 15:
            SetViewExtentsType(buf.ReadInt());
            break;
        case 16:
            SetViewDimension(buf.ReadInt());
            break;
        case 17:
            SetViewKeyframes(buf.ReadIntVector());
            break;
        case 18:
            SetCameraViewMode(buf.ReadBool());
            break;
        case 19:
            SetUsingScalableRendering(buf.ReadBool());
            break;
        case 20:
            SetLastRenderMin(buf.ReadFloat());
            break;
        case 21:
            SetLastRenderAvg(buf.ReadFloat());
            break;
        case 22:
            SetLastRenderMax(buf.ReadFloat());
            break;
        case 23:
            SetNumPrimitives(buf.ReadInt());
            break;
        case 24:
            SetExtents(buf.ReadDoubleArray());
            break;
        case 25:
            SetWindowSize(buf.ReadIntArray());
            break;
        case 26:
            SetWinMode(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("activeSource", activeSource, indent) + "\n";
        str = str + intToString("activeTimeSlider", activeTimeSlider, indent) + "\n";
        str = str + stringVectorToString("timeSliders", timeSliders, indent) + "\n";
        str = str + intVectorToString("timeSliderCurrentStates", timeSliderCurrentStates, indent) + "\n";
        str = str + intToString("animationMode", animationMode, indent) + "\n";
        str = str + intToString("interactionMode", interactionMode, indent) + "\n";
        str = str + intToString("toolUpdateMode", toolUpdateMode, indent) + "\n";
        str = str + boolToString("boundingBoxNavigate", boundingBoxNavigate, indent) + "\n";
        str = str + boolToString("spin", spin, indent) + "\n";
        str = str + boolToString("fullFrame", fullFrame, indent) + "\n";
        str = str + boolToString("perspective", perspective, indent) + "\n";
        str = str + boolToString("maintainView", maintainView, indent) + "\n";
        str = str + boolToString("lockView", lockView, indent) + "\n";
        str = str + boolToString("lockTools", lockTools, indent) + "\n";
        str = str + boolToString("lockTime", lockTime, indent) + "\n";
        str = str + intToString("viewExtentsType", viewExtentsType, indent) + "\n";
        str = str + intToString("viewDimension", viewDimension, indent) + "\n";
        str = str + intVectorToString("viewKeyframes", viewKeyframes, indent) + "\n";
        str = str + boolToString("cameraViewMode", cameraViewMode, indent) + "\n";
        str = str + boolToString("usingScalableRendering", usingScalableRendering, indent) + "\n";
        str = str + floatToString("lastRenderMin", lastRenderMin, indent) + "\n";
        str = str + floatToString("lastRenderAvg", lastRenderAvg, indent) + "\n";
        str = str + floatToString("lastRenderMax", lastRenderMax, indent) + "\n";
        str = str + intToString("numPrimitives", numPrimitives, indent) + "\n";
        str = str + doubleArrayToString("extents", extents, indent) + "\n";
        str = str + intArrayToString("windowSize", windowSize, indent) + "\n";
        str = str + intToString("winMode", winMode, indent) + "\n";
        return str;
    }


    // Attributes
    private String   activeSource;
    private int      activeTimeSlider;
    private Vector   timeSliders; // vector of String objects
    private Vector   timeSliderCurrentStates; // vector of Integer objects
    private int      animationMode;
    private int      interactionMode;
    private int      toolUpdateMode;
    private boolean  boundingBoxNavigate;
    private boolean  spin;
    private boolean  fullFrame;
    private boolean  perspective;
    private boolean  maintainView;
    private boolean  lockView;
    private boolean  lockTools;
    private boolean  lockTime;
    private int      viewExtentsType;
    private int      viewDimension;
    private Vector   viewKeyframes; // vector of Integer objects
    private boolean  cameraViewMode;
    private boolean  usingScalableRendering;
    private float    lastRenderMin;
    private float    lastRenderAvg;
    private float    lastRenderMax;
    private int      numPrimitives;
    private double[] extents;
    private int[]    windowSize;
    private int      winMode;
}

