/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef GLOBALATTRIBUTES_H
#define GLOBALATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: GlobalAttributes
//
// Purpose:
//    This class contains attributes associated with the main window.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Mon Feb 4 15:22:16 PST 2008
//
// Modifications:
//   
// ****************************************************************************

class STATE_API GlobalAttributes : public AttributeSubject
{
public:
    GlobalAttributes();
    GlobalAttributes(const GlobalAttributes &obj);
    virtual ~GlobalAttributes();

    virtual GlobalAttributes& operator = (const GlobalAttributes &obj);
    virtual bool operator == (const GlobalAttributes &obj) const;
    virtual bool operator != (const GlobalAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectSources();
    void SelectWindows();

    // Property setting methods
    void SetSources(const stringVector &sources_);
    void SetWindows(const intVector &windows_);
    void SetActiveWindow(int activeWindow_);
    void SetIconifiedFlag(bool iconifiedFlag_);
    void SetAutoUpdateFlag(bool autoUpdateFlag_);
    void SetReplacePlots(bool replacePlots_);
    void SetApplyOperator(bool applyOperator_);
    void SetExecuting(bool executing_);
    void SetWindowLayout(int windowLayout_);
    void SetMakeDefaultConfirm(bool makeDefaultConfirm_);
    void SetCloneWindowOnFirstRef(bool cloneWindowOnFirstRef_);
    void SetMaintainView(bool maintainView_);
    void SetMaintainData(bool maintainData_);
    void SetAutomaticallyAddOperator(bool automaticallyAddOperator_);
    void SetTryHarderCyclesTimes(bool tryHarderCyclesTimes_);
    void SetTreatAllDBsAsTimeVarying(bool treatAllDBsAsTimeVarying_);
    void SetCreateMeshQualityExpressions(bool createMeshQualityExpressions_);
    void SetCreateTimeDerivativeExpressions(bool createTimeDerivativeExpressions_);
    void SetCreateVectorMagnitudeExpressions(bool createVectorMagnitudeExpressions_);
    void SetNewPlotsInheritSILRestriction(bool newPlotsInheritSILRestriction_);
    void SetUserDirForSessionFiles(bool userDirForSessionFiles_);
    void SetSaveCrashRecoveryFile(bool saveCrashRecoveryFile_);
    void SetApplySelection(bool applySelection_);

    // Property getting methods
    const stringVector &GetSources() const;
          stringVector &GetSources();
    const intVector    &GetWindows() const;
          intVector    &GetWindows();
    int                GetActiveWindow() const;
    bool               GetIconifiedFlag() const;
    bool               GetAutoUpdateFlag() const;
    bool               GetReplacePlots() const;
    bool               GetApplyOperator() const;
    bool               GetExecuting() const;
    int                GetWindowLayout() const;
    bool               GetMakeDefaultConfirm() const;
    bool               GetCloneWindowOnFirstRef() const;
    bool               GetMaintainView() const;
    bool               GetMaintainData() const;
    bool               GetAutomaticallyAddOperator() const;
    bool               GetTryHarderCyclesTimes() const;
    bool               GetTreatAllDBsAsTimeVarying() const;
    bool               GetCreateMeshQualityExpressions() const;
    bool               GetCreateTimeDerivativeExpressions() const;
    bool               GetCreateVectorMagnitudeExpressions() const;
    bool               GetNewPlotsInheritSILRestriction() const;
    bool               GetUserDirForSessionFiles() const;
    bool               GetSaveCrashRecoveryFile() const;
    bool               GetApplySelection() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_sources = 0,
        ID_windows,
        ID_activeWindow,
        ID_iconifiedFlag,
        ID_autoUpdateFlag,
        ID_replacePlots,
        ID_applyOperator,
        ID_executing,
        ID_windowLayout,
        ID_makeDefaultConfirm,
        ID_cloneWindowOnFirstRef,
        ID_maintainView,
        ID_maintainData,
        ID_automaticallyAddOperator,
        ID_tryHarderCyclesTimes,
        ID_treatAllDBsAsTimeVarying,
        ID_createMeshQualityExpressions,
        ID_createTimeDerivativeExpressions,
        ID_createVectorMagnitudeExpressions,
        ID_newPlotsInheritSILRestriction,
        ID_userDirForSessionFiles,
        ID_saveCrashRecoveryFile,
        ID_applySelection
    };

private:
    stringVector sources;
    intVector    windows;
    int          activeWindow;
    bool         iconifiedFlag;
    bool         autoUpdateFlag;
    bool         replacePlots;
    bool         applyOperator;
    bool         executing;
    int          windowLayout;
    bool         makeDefaultConfirm;
    bool         cloneWindowOnFirstRef;
    bool         maintainView;
    bool         maintainData;
    bool         automaticallyAddOperator;
    bool         tryHarderCyclesTimes;
    bool         treatAllDBsAsTimeVarying;
    bool         createMeshQualityExpressions;
    bool         createTimeDerivativeExpressions;
    bool         createVectorMagnitudeExpressions;
    bool         newPlotsInheritSILRestriction;
    bool         userDirForSessionFiles;
    bool         saveCrashRecoveryFile;
    bool         applySelection;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
