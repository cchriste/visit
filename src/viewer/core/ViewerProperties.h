/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
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

#ifndef VIEWERPROPERTIES_H
#define VIEWERPROPERTIES_H
#include <viewercore_exports.h>
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: ViewerProperties
//
// Purpose:
//    Contain properties about the viewer
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class VIEWERCORE_API ViewerProperties : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    ViewerProperties();
    ViewerProperties(const ViewerProperties &obj);
protected:
    // These constructors are for objects derived from this class
    ViewerProperties(private_tmfs_t tmfs);
    ViewerProperties(const ViewerProperties &obj, private_tmfs_t tmfs);
public:
    virtual ~ViewerProperties();

    virtual ViewerProperties& operator = (const ViewerProperties &obj);
    virtual bool operator == (const ViewerProperties &obj) const;
    virtual bool operator != (const ViewerProperties &obj) const;
private:
    void Init();
    void Copy(const ViewerProperties &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectWindowBorders();
    void SelectWindowShift();
    void SelectWindowPreShift();
    void SelectWindowGeometry();
    void SelectConfigurationFileName();
    void SelectApplicationLocale();
    void SelectEngineParallelArguments();
    void SelectAssumedFormats();
    void SelectFallbackFormats();

    // Property setting methods
    void SetMasterProcess(bool MasterProcess_);
    void SetInSitu(bool InSitu_);
    void SetNowin(bool nowin_);
    void SetWindowBorders(const std::string &windowBorders_);
    void SetWindowShift(const std::string &windowShift_);
    void SetWindowPreShift(const std::string &windowPreShift_);
    void SetWindowGeometry(const std::string &windowGeometry_);
    void SetWindowFullScreen(int windowFullScreen_);
    void SetWindowSmall(bool windowSmall_);
    void SetNoConfig(bool noConfig_);
    void SetDefaultStereoToOn(bool defaultStereoToOn_);
    void SetUseWindowMetrics(bool useWindowMetrics_);
    void SetConfigurationFileName(const std::string &configurationFileName_);
    void SetApplicationLocale(const std::string &applicationLocale_);
    void SetLaunchedByClient(bool launchedByClient_);
    void SetSuppressMessages(bool suppressMessages_);
    void SetNumEngineRestarts(int numEngineRestarts_);
    void SetEngineParallelArguments(const stringVector &engineParallelArguments_);
    void SetAssumedFormats(const stringVector &assumedFormats_);
    void SetFallbackFormats(const stringVector &fallbackFormats_);
    void SetDebugLevel(int debugLevel_);
    void SetBufferDebug(bool bufferDebug_);
    void SetForceSSHTunneling(bool forceSSHTunneling_);
    void SetInExecute(bool inExecute_);
    void SetInLaunch(bool inLaunch_);
    void SetDecorateDebug(bool decorateDebug_);

    // Property getting methods
    bool               GetMasterProcess() const;
    bool               GetInSitu() const;
    bool               GetNowin() const;
    const std::string  &GetWindowBorders() const;
          std::string  &GetWindowBorders();
    const std::string  &GetWindowShift() const;
          std::string  &GetWindowShift();
    const std::string  &GetWindowPreShift() const;
          std::string  &GetWindowPreShift();
    const std::string  &GetWindowGeometry() const;
          std::string  &GetWindowGeometry();
    int                GetWindowFullScreen() const;
    bool               GetWindowSmall() const;
    bool               GetNoConfig() const;
    bool               GetDefaultStereoToOn() const;
    bool               GetUseWindowMetrics() const;
    const std::string  &GetConfigurationFileName() const;
          std::string  &GetConfigurationFileName();
    const std::string  &GetApplicationLocale() const;
          std::string  &GetApplicationLocale();
    bool               GetLaunchedByClient() const;
    bool               GetSuppressMessages() const;
    int                GetNumEngineRestarts() const;
    const stringVector &GetEngineParallelArguments() const;
          stringVector &GetEngineParallelArguments();
    const stringVector &GetAssumedFormats() const;
          stringVector &GetAssumedFormats();
    const stringVector &GetFallbackFormats() const;
          stringVector &GetFallbackFormats();
    int                GetDebugLevel() const;
    bool               GetBufferDebug() const;
    bool               GetForceSSHTunneling() const;
    bool               GetInExecute() const;
    bool               GetInLaunch() const;
    bool               GetDecorateDebug() const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_MasterProcess = 0,
        ID_InSitu,
        ID_nowin,
        ID_windowBorders,
        ID_windowShift,
        ID_windowPreShift,
        ID_windowGeometry,
        ID_windowFullScreen,
        ID_windowSmall,
        ID_noConfig,
        ID_defaultStereoToOn,
        ID_useWindowMetrics,
        ID_configurationFileName,
        ID_applicationLocale,
        ID_launchedByClient,
        ID_suppressMessages,
        ID_numEngineRestarts,
        ID_engineParallelArguments,
        ID_assumedFormats,
        ID_fallbackFormats,
        ID_debugLevel,
        ID_bufferDebug,
        ID_forceSSHTunneling,
        ID_inExecute,
        ID_inLaunch,
        ID_decorateDebug,
        ID__LAST
    };

private:
    bool         MasterProcess;
    bool         InSitu;
    bool         nowin;
    std::string  windowBorders;
    std::string  windowShift;
    std::string  windowPreShift;
    std::string  windowGeometry;
    int          windowFullScreen;
    bool         windowSmall;
    bool         noConfig;
    bool         defaultStereoToOn;
    bool         useWindowMetrics;
    std::string  configurationFileName;
    std::string  applicationLocale;
    bool         launchedByClient;
    bool         suppressMessages;
    int          numEngineRestarts;
    stringVector engineParallelArguments;
    stringVector assumedFormats;
    stringVector fallbackFormats;
    int          debugLevel;
    bool         bufferDebug;
    bool         forceSSHTunneling;
    bool         inExecute;
    bool         inLaunch;
    bool         decorateDebug;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define VIEWERPROPERTIES_TMFS "bbbssssibbbbssbbis*s*s*ibbbbb"

#endif
