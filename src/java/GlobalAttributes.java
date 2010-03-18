// ***************************************************************************
//
// Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400124
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
// Class: GlobalAttributes
//
// Purpose:
//    This class contains attributes associated with the main window.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class GlobalAttributes extends AttributeSubject
{
    private static int numAdditionalAttributes = 23;

    public GlobalAttributes()
    {
        super(numAdditionalAttributes);

        sources = new Vector();
        windows = new Vector();
        activeWindow = 1;
        iconifiedFlag = false;
        autoUpdateFlag = false;
        replacePlots = false;
        applyOperator = true;
        executing = false;
        windowLayout = 1;
        makeDefaultConfirm = true;
        cloneWindowOnFirstRef = false;
        automaticallyAddOperator = false;
        tryHarderCyclesTimes = false;
        treatAllDBsAsTimeVarying = false;
        createMeshQualityExpressions = true;
        createTimeDerivativeExpressions = true;
        createVectorMagnitudeExpressions = true;
        newPlotsInheritSILRestriction = true;
        userDirForSessionFiles = false;
        saveCrashRecoveryFile = true;
        applySelection = true;
        ignoreExtentsFromDbs = false;
        expandNewPlots = false;
    }

    public GlobalAttributes(int nMoreFields)
    {
        super(numAdditionalAttributes + nMoreFields);

        sources = new Vector();
        windows = new Vector();
        activeWindow = 1;
        iconifiedFlag = false;
        autoUpdateFlag = false;
        replacePlots = false;
        applyOperator = true;
        executing = false;
        windowLayout = 1;
        makeDefaultConfirm = true;
        cloneWindowOnFirstRef = false;
        automaticallyAddOperator = false;
        tryHarderCyclesTimes = false;
        treatAllDBsAsTimeVarying = false;
        createMeshQualityExpressions = true;
        createTimeDerivativeExpressions = true;
        createVectorMagnitudeExpressions = true;
        newPlotsInheritSILRestriction = true;
        userDirForSessionFiles = false;
        saveCrashRecoveryFile = true;
        applySelection = true;
        ignoreExtentsFromDbs = false;
        expandNewPlots = false;
    }

    public GlobalAttributes(GlobalAttributes obj)
    {
        super(numAdditionalAttributes);

        int i;

        sources = new Vector(obj.sources.size());
        for(i = 0; i < obj.sources.size(); ++i)
            sources.addElement(new String((String)obj.sources.elementAt(i)));

        windows = new Vector();
        for(i = 0; i < obj.windows.size(); ++i)
        {
            Integer iv = (Integer)obj.windows.elementAt(i);
            windows.addElement(new Integer(iv.intValue()));
        }
        activeWindow = obj.activeWindow;
        iconifiedFlag = obj.iconifiedFlag;
        autoUpdateFlag = obj.autoUpdateFlag;
        replacePlots = obj.replacePlots;
        applyOperator = obj.applyOperator;
        executing = obj.executing;
        windowLayout = obj.windowLayout;
        makeDefaultConfirm = obj.makeDefaultConfirm;
        cloneWindowOnFirstRef = obj.cloneWindowOnFirstRef;
        automaticallyAddOperator = obj.automaticallyAddOperator;
        tryHarderCyclesTimes = obj.tryHarderCyclesTimes;
        treatAllDBsAsTimeVarying = obj.treatAllDBsAsTimeVarying;
        createMeshQualityExpressions = obj.createMeshQualityExpressions;
        createTimeDerivativeExpressions = obj.createTimeDerivativeExpressions;
        createVectorMagnitudeExpressions = obj.createVectorMagnitudeExpressions;
        newPlotsInheritSILRestriction = obj.newPlotsInheritSILRestriction;
        userDirForSessionFiles = obj.userDirForSessionFiles;
        saveCrashRecoveryFile = obj.saveCrashRecoveryFile;
        applySelection = obj.applySelection;
        ignoreExtentsFromDbs = obj.ignoreExtentsFromDbs;
        expandNewPlots = obj.expandNewPlots;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return numAdditionalAttributes;
    }

    public boolean equals(GlobalAttributes obj)
    {
        int i;

        // Compare the elements in the sources vector.
        boolean sources_equal = (obj.sources.size() == sources.size());
        for(i = 0; (i < sources.size()) && sources_equal; ++i)
        {
            // Make references to String from Object.
            String sources1 = (String)sources.elementAt(i);
            String sources2 = (String)obj.sources.elementAt(i);
            sources_equal = sources1.equals(sources2);
        }
        // Compare the elements in the windows vector.
        boolean windows_equal = (obj.windows.size() == windows.size());
        for(i = 0; (i < windows.size()) && windows_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer windows1 = (Integer)windows.elementAt(i);
            Integer windows2 = (Integer)obj.windows.elementAt(i);
            windows_equal = windows1.equals(windows2);
        }
        // Create the return value
        return (sources_equal &&
                windows_equal &&
                (activeWindow == obj.activeWindow) &&
                (iconifiedFlag == obj.iconifiedFlag) &&
                (autoUpdateFlag == obj.autoUpdateFlag) &&
                (replacePlots == obj.replacePlots) &&
                (applyOperator == obj.applyOperator) &&
                (executing == obj.executing) &&
                (windowLayout == obj.windowLayout) &&
                (makeDefaultConfirm == obj.makeDefaultConfirm) &&
                (cloneWindowOnFirstRef == obj.cloneWindowOnFirstRef) &&
                (automaticallyAddOperator == obj.automaticallyAddOperator) &&
                (tryHarderCyclesTimes == obj.tryHarderCyclesTimes) &&
                (treatAllDBsAsTimeVarying == obj.treatAllDBsAsTimeVarying) &&
                (createMeshQualityExpressions == obj.createMeshQualityExpressions) &&
                (createTimeDerivativeExpressions == obj.createTimeDerivativeExpressions) &&
                (createVectorMagnitudeExpressions == obj.createVectorMagnitudeExpressions) &&
                (newPlotsInheritSILRestriction == obj.newPlotsInheritSILRestriction) &&
                (userDirForSessionFiles == obj.userDirForSessionFiles) &&
                (saveCrashRecoveryFile == obj.saveCrashRecoveryFile) &&
                (applySelection == obj.applySelection) &&
                (ignoreExtentsFromDbs == obj.ignoreExtentsFromDbs) &&
                (expandNewPlots == obj.expandNewPlots));
    }

    // Property setting methods
    public void SetSources(Vector sources_)
    {
        sources = sources_;
        Select(0);
    }

    public void SetWindows(Vector windows_)
    {
        windows = windows_;
        Select(1);
    }

    public void SetActiveWindow(int activeWindow_)
    {
        activeWindow = activeWindow_;
        Select(2);
    }

    public void SetIconifiedFlag(boolean iconifiedFlag_)
    {
        iconifiedFlag = iconifiedFlag_;
        Select(3);
    }

    public void SetAutoUpdateFlag(boolean autoUpdateFlag_)
    {
        autoUpdateFlag = autoUpdateFlag_;
        Select(4);
    }

    public void SetReplacePlots(boolean replacePlots_)
    {
        replacePlots = replacePlots_;
        Select(5);
    }

    public void SetApplyOperator(boolean applyOperator_)
    {
        applyOperator = applyOperator_;
        Select(6);
    }

    public void SetExecuting(boolean executing_)
    {
        executing = executing_;
        Select(7);
    }

    public void SetWindowLayout(int windowLayout_)
    {
        windowLayout = windowLayout_;
        Select(8);
    }

    public void SetMakeDefaultConfirm(boolean makeDefaultConfirm_)
    {
        makeDefaultConfirm = makeDefaultConfirm_;
        Select(9);
    }

    public void SetCloneWindowOnFirstRef(boolean cloneWindowOnFirstRef_)
    {
        cloneWindowOnFirstRef = cloneWindowOnFirstRef_;
        Select(10);
    }

    public void SetAutomaticallyAddOperator(boolean automaticallyAddOperator_)
    {
        automaticallyAddOperator = automaticallyAddOperator_;
        Select(11);
    }

    public void SetTryHarderCyclesTimes(boolean tryHarderCyclesTimes_)
    {
        tryHarderCyclesTimes = tryHarderCyclesTimes_;
        Select(12);
    }

    public void SetTreatAllDBsAsTimeVarying(boolean treatAllDBsAsTimeVarying_)
    {
        treatAllDBsAsTimeVarying = treatAllDBsAsTimeVarying_;
        Select(13);
    }

    public void SetCreateMeshQualityExpressions(boolean createMeshQualityExpressions_)
    {
        createMeshQualityExpressions = createMeshQualityExpressions_;
        Select(14);
    }

    public void SetCreateTimeDerivativeExpressions(boolean createTimeDerivativeExpressions_)
    {
        createTimeDerivativeExpressions = createTimeDerivativeExpressions_;
        Select(15);
    }

    public void SetCreateVectorMagnitudeExpressions(boolean createVectorMagnitudeExpressions_)
    {
        createVectorMagnitudeExpressions = createVectorMagnitudeExpressions_;
        Select(16);
    }

    public void SetNewPlotsInheritSILRestriction(boolean newPlotsInheritSILRestriction_)
    {
        newPlotsInheritSILRestriction = newPlotsInheritSILRestriction_;
        Select(17);
    }

    public void SetUserDirForSessionFiles(boolean userDirForSessionFiles_)
    {
        userDirForSessionFiles = userDirForSessionFiles_;
        Select(18);
    }

    public void SetSaveCrashRecoveryFile(boolean saveCrashRecoveryFile_)
    {
        saveCrashRecoveryFile = saveCrashRecoveryFile_;
        Select(19);
    }

    public void SetApplySelection(boolean applySelection_)
    {
        applySelection = applySelection_;
        Select(20);
    }

    public void SetIgnoreExtentsFromDbs(boolean ignoreExtentsFromDbs_)
    {
        ignoreExtentsFromDbs = ignoreExtentsFromDbs_;
        Select(21);
    }

    public void SetExpandNewPlots(boolean expandNewPlots_)
    {
        expandNewPlots = expandNewPlots_;
        Select(22);
    }

    // Property getting methods
    public Vector  GetSources() { return sources; }
    public Vector  GetWindows() { return windows; }
    public int     GetActiveWindow() { return activeWindow; }
    public boolean GetIconifiedFlag() { return iconifiedFlag; }
    public boolean GetAutoUpdateFlag() { return autoUpdateFlag; }
    public boolean GetReplacePlots() { return replacePlots; }
    public boolean GetApplyOperator() { return applyOperator; }
    public boolean GetExecuting() { return executing; }
    public int     GetWindowLayout() { return windowLayout; }
    public boolean GetMakeDefaultConfirm() { return makeDefaultConfirm; }
    public boolean GetCloneWindowOnFirstRef() { return cloneWindowOnFirstRef; }
    public boolean GetAutomaticallyAddOperator() { return automaticallyAddOperator; }
    public boolean GetTryHarderCyclesTimes() { return tryHarderCyclesTimes; }
    public boolean GetTreatAllDBsAsTimeVarying() { return treatAllDBsAsTimeVarying; }
    public boolean GetCreateMeshQualityExpressions() { return createMeshQualityExpressions; }
    public boolean GetCreateTimeDerivativeExpressions() { return createTimeDerivativeExpressions; }
    public boolean GetCreateVectorMagnitudeExpressions() { return createVectorMagnitudeExpressions; }
    public boolean GetNewPlotsInheritSILRestriction() { return newPlotsInheritSILRestriction; }
    public boolean GetUserDirForSessionFiles() { return userDirForSessionFiles; }
    public boolean GetSaveCrashRecoveryFile() { return saveCrashRecoveryFile; }
    public boolean GetApplySelection() { return applySelection; }
    public boolean GetIgnoreExtentsFromDbs() { return ignoreExtentsFromDbs; }
    public boolean GetExpandNewPlots() { return expandNewPlots; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteStringVector(sources);
        if(WriteSelect(1, buf))
            buf.WriteIntVector(windows);
        if(WriteSelect(2, buf))
            buf.WriteInt(activeWindow);
        if(WriteSelect(3, buf))
            buf.WriteBool(iconifiedFlag);
        if(WriteSelect(4, buf))
            buf.WriteBool(autoUpdateFlag);
        if(WriteSelect(5, buf))
            buf.WriteBool(replacePlots);
        if(WriteSelect(6, buf))
            buf.WriteBool(applyOperator);
        if(WriteSelect(7, buf))
            buf.WriteBool(executing);
        if(WriteSelect(8, buf))
            buf.WriteInt(windowLayout);
        if(WriteSelect(9, buf))
            buf.WriteBool(makeDefaultConfirm);
        if(WriteSelect(10, buf))
            buf.WriteBool(cloneWindowOnFirstRef);
        if(WriteSelect(11, buf))
            buf.WriteBool(automaticallyAddOperator);
        if(WriteSelect(12, buf))
            buf.WriteBool(tryHarderCyclesTimes);
        if(WriteSelect(13, buf))
            buf.WriteBool(treatAllDBsAsTimeVarying);
        if(WriteSelect(14, buf))
            buf.WriteBool(createMeshQualityExpressions);
        if(WriteSelect(15, buf))
            buf.WriteBool(createTimeDerivativeExpressions);
        if(WriteSelect(16, buf))
            buf.WriteBool(createVectorMagnitudeExpressions);
        if(WriteSelect(17, buf))
            buf.WriteBool(newPlotsInheritSILRestriction);
        if(WriteSelect(18, buf))
            buf.WriteBool(userDirForSessionFiles);
        if(WriteSelect(19, buf))
            buf.WriteBool(saveCrashRecoveryFile);
        if(WriteSelect(20, buf))
            buf.WriteBool(applySelection);
        if(WriteSelect(21, buf))
            buf.WriteBool(ignoreExtentsFromDbs);
        if(WriteSelect(22, buf))
            buf.WriteBool(expandNewPlots);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetSources(buf.ReadStringVector());
            break;
        case 1:
            SetWindows(buf.ReadIntVector());
            break;
        case 2:
            SetActiveWindow(buf.ReadInt());
            break;
        case 3:
            SetIconifiedFlag(buf.ReadBool());
            break;
        case 4:
            SetAutoUpdateFlag(buf.ReadBool());
            break;
        case 5:
            SetReplacePlots(buf.ReadBool());
            break;
        case 6:
            SetApplyOperator(buf.ReadBool());
            break;
        case 7:
            SetExecuting(buf.ReadBool());
            break;
        case 8:
            SetWindowLayout(buf.ReadInt());
            break;
        case 9:
            SetMakeDefaultConfirm(buf.ReadBool());
            break;
        case 10:
            SetCloneWindowOnFirstRef(buf.ReadBool());
            break;
        case 11:
            SetAutomaticallyAddOperator(buf.ReadBool());
            break;
        case 12:
            SetTryHarderCyclesTimes(buf.ReadBool());
            break;
        case 13:
            SetTreatAllDBsAsTimeVarying(buf.ReadBool());
            break;
        case 14:
            SetCreateMeshQualityExpressions(buf.ReadBool());
            break;
        case 15:
            SetCreateTimeDerivativeExpressions(buf.ReadBool());
            break;
        case 16:
            SetCreateVectorMagnitudeExpressions(buf.ReadBool());
            break;
        case 17:
            SetNewPlotsInheritSILRestriction(buf.ReadBool());
            break;
        case 18:
            SetUserDirForSessionFiles(buf.ReadBool());
            break;
        case 19:
            SetSaveCrashRecoveryFile(buf.ReadBool());
            break;
        case 20:
            SetApplySelection(buf.ReadBool());
            break;
        case 21:
            SetIgnoreExtentsFromDbs(buf.ReadBool());
            break;
        case 22:
            SetExpandNewPlots(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringVectorToString("sources", sources, indent) + "\n";
        str = str + intVectorToString("windows", windows, indent) + "\n";
        str = str + intToString("activeWindow", activeWindow, indent) + "\n";
        str = str + boolToString("iconifiedFlag", iconifiedFlag, indent) + "\n";
        str = str + boolToString("autoUpdateFlag", autoUpdateFlag, indent) + "\n";
        str = str + boolToString("replacePlots", replacePlots, indent) + "\n";
        str = str + boolToString("applyOperator", applyOperator, indent) + "\n";
        str = str + boolToString("executing", executing, indent) + "\n";
        str = str + intToString("windowLayout", windowLayout, indent) + "\n";
        str = str + boolToString("makeDefaultConfirm", makeDefaultConfirm, indent) + "\n";
        str = str + boolToString("cloneWindowOnFirstRef", cloneWindowOnFirstRef, indent) + "\n";
        str = str + boolToString("automaticallyAddOperator", automaticallyAddOperator, indent) + "\n";
        str = str + boolToString("tryHarderCyclesTimes", tryHarderCyclesTimes, indent) + "\n";
        str = str + boolToString("treatAllDBsAsTimeVarying", treatAllDBsAsTimeVarying, indent) + "\n";
        str = str + boolToString("createMeshQualityExpressions", createMeshQualityExpressions, indent) + "\n";
        str = str + boolToString("createTimeDerivativeExpressions", createTimeDerivativeExpressions, indent) + "\n";
        str = str + boolToString("createVectorMagnitudeExpressions", createVectorMagnitudeExpressions, indent) + "\n";
        str = str + boolToString("newPlotsInheritSILRestriction", newPlotsInheritSILRestriction, indent) + "\n";
        str = str + boolToString("userDirForSessionFiles", userDirForSessionFiles, indent) + "\n";
        str = str + boolToString("saveCrashRecoveryFile", saveCrashRecoveryFile, indent) + "\n";
        str = str + boolToString("applySelection", applySelection, indent) + "\n";
        str = str + boolToString("ignoreExtentsFromDbs", ignoreExtentsFromDbs, indent) + "\n";
        str = str + boolToString("expandNewPlots", expandNewPlots, indent) + "\n";
        return str;
    }


    // Attributes
    private Vector  sources; // vector of String objects
    private Vector  windows; // vector of Integer objects
    private int     activeWindow;
    private boolean iconifiedFlag;
    private boolean autoUpdateFlag;
    private boolean replacePlots;
    private boolean applyOperator;
    private boolean executing;
    private int     windowLayout;
    private boolean makeDefaultConfirm;
    private boolean cloneWindowOnFirstRef;
    private boolean automaticallyAddOperator;
    private boolean tryHarderCyclesTimes;
    private boolean treatAllDBsAsTimeVarying;
    private boolean createMeshQualityExpressions;
    private boolean createTimeDerivativeExpressions;
    private boolean createVectorMagnitudeExpressions;
    private boolean newPlotsInheritSILRestriction;
    private boolean userDirForSessionFiles;
    private boolean saveCrashRecoveryFile;
    private boolean applySelection;
    private boolean ignoreExtentsFromDbs;
    private boolean expandNewPlots;
}

