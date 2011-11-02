/*****************************************************************************
*
* Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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

#ifndef LINESAMPLERATTRIBUTES_H
#define LINESAMPLERATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: LineSamplerAttributes
//
// Purpose:
//    This class contains attributes for the line sampler operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class LineSamplerAttributes : public AttributeSubject
{
public:
    enum MeshGeometry
    {
        Cartesian,
        Cylindrical,
        Toroidal
    };
    enum ArrayConfiguration
    {
        Geometry,
        Manual
    };
    enum Boundary
    {
        Data,
        Wall
    };
    enum ChannelProjection
    {
        Divergent,
        Parallel,
        Grid
    };
    enum ChannelLayoutType
    {
        ChannelAbsolute,
        ChannelRelative
    };
    enum ArrayAxis
    {
        R,
        Z
    };
    enum ViewGeometry
    {
        Points,
        Lines,
        Surfaces
    };
    enum DisplayTime
    {
        Step,
        Time,
        Cycle
    };
    enum ChannelGeometry
    {
        Point,
        Line,
        Cylinder,
        Cone
    };
    enum ViewDimension
    {
        One,
        Two,
        Three
    };
    enum ChannelProfile
    {
        TopHat,
        Gaussian
    };
    enum ChannelIntegration
    {
        NoChannelIntegration,
        IntegrateAlongChannel
    };
    enum ToroidalIntegration
    {
        NoToroidalIntegration,
        ToroidalTimeSample,
        IntegrateToroidally
    };
    enum ToroidalAngleSampling
    {
        ToroidalAngleAbsoluteSampling,
        ToroidalAngleRelativeSampling
    };
    enum TimeSampling
    {
        CurrentTimeStep,
        MultipleTimeSteps
    };

    // These constructors are for objects of this class
    LineSamplerAttributes();
    LineSamplerAttributes(const LineSamplerAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    LineSamplerAttributes(private_tmfs_t tmfs);
    LineSamplerAttributes(const LineSamplerAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~LineSamplerAttributes();

    virtual LineSamplerAttributes& operator = (const LineSamplerAttributes &obj);
    virtual bool operator == (const LineSamplerAttributes &obj) const;
    virtual bool operator != (const LineSamplerAttributes &obj) const;
private:
    void Init();
    void Copy(const LineSamplerAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectArrayOrigin();
    void SelectChannelList();
    void SelectWallList();

    // Property setting methods
    void SetMeshGeometry(MeshGeometry meshGeometry_);
    void SetArrayConfiguration(ArrayConfiguration arrayConfiguration_);
    void SetBoundary(Boundary boundary_);
    void SetNArrays(int nArrays_);
    void SetToroidalArrayAngle(double toroidalArrayAngle_);
    void SetNChannels(int nChannels_);
    void SetChannelProjection(ChannelProjection channelProjection_);
    void SetChannelLayoutType(ChannelLayoutType channelLayoutType_);
    void SetChannelOffset(double channelOffset_);
    void SetChannelAngle(double channelAngle_);
    void SetNRows(int nRows_);
    void SetRowOffset(double rowOffset_);
    void SetArrayOrigin(const double *arrayOrigin_);
    void SetArrayAxis(ArrayAxis arrayAxis_);
    void SetPoloidalAngleStart(double poloidalAngleStart_);
    void SetPoloidalAngleStop(double poloidalAngleStop_);
    void SetPoloialAngle(double poloialAngle_);
    void SetPoloialRTilt(double poloialRTilt_);
    void SetPoloialZTilt(double poloialZTilt_);
    void SetToroidalAngle(double toroidalAngle_);
    void SetFlipToroidalAngle(bool flipToroidalAngle_);
    void SetViewGeometry(ViewGeometry viewGeometry_);
    void SetViewDimension(ViewDimension viewDimension_);
    void SetHeightPlotScale(double heightPlotScale_);
    void SetChannelPlotOffset(double channelPlotOffset_);
    void SetArrayPlotOffset(double arrayPlotOffset_);
    void SetDisplayTime(DisplayTime displayTime_);
    void SetChannelGeometry(ChannelGeometry channelGeometry_);
    void SetRadius(double radius_);
    void SetDivergence(double divergence_);
    void SetChannelProfile(ChannelProfile channelProfile_);
    void SetStandardDeviation(double standardDeviation_);
    void SetSampleDistance(double sampleDistance_);
    void SetSampleVolume(double sampleVolume_);
    void SetSampleArc(double sampleArc_);
    void SetChannelIntegration(ChannelIntegration channelIntegration_);
    void SetToroidalIntegration(ToroidalIntegration toroidalIntegration_);
    void SetToroidalAngleSampling(ToroidalAngleSampling toroidalAngleSampling_);
    void SetToroidalAngleStart(double toroidalAngleStart_);
    void SetToroidalAngleStop(double toroidalAngleStop_);
    void SetToroidalAngleStride(double toroidalAngleStride_);
    void SetTimeSampling(TimeSampling timeSampling_);
    void SetTimeStepStart(int timeStepStart_);
    void SetTimeStepStop(int timeStepStop_);
    void SetTimeStepStride(int timeStepStride_);
    void SetChannelList(const doubleVector &channelList_);
    void SetWallList(const doubleVector &wallList_);
    void SetNChannelListArrays(int nChannelListArrays_);
    void SetChannelListToroidalArrayAngle(double channelListToroidalArrayAngle_);
    void SetChannelListToroidalAngle(double channelListToroidalAngle_);

    // Property getting methods
    MeshGeometry       GetMeshGeometry() const;
    ArrayConfiguration GetArrayConfiguration() const;
    Boundary           GetBoundary() const;
    int                GetNArrays() const;
    double             GetToroidalArrayAngle() const;
    int                GetNChannels() const;
    ChannelProjection  GetChannelProjection() const;
    ChannelLayoutType  GetChannelLayoutType() const;
    double             GetChannelOffset() const;
    double             GetChannelAngle() const;
    int                GetNRows() const;
    double             GetRowOffset() const;
    const double       *GetArrayOrigin() const;
          double       *GetArrayOrigin();
    ArrayAxis          GetArrayAxis() const;
    double             GetPoloidalAngleStart() const;
    double             GetPoloidalAngleStop() const;
    double             GetPoloialAngle() const;
    double             GetPoloialRTilt() const;
    double             GetPoloialZTilt() const;
    double             GetToroidalAngle() const;
    bool               GetFlipToroidalAngle() const;
    ViewGeometry       GetViewGeometry() const;
    ViewDimension      GetViewDimension() const;
    double             GetHeightPlotScale() const;
    double             GetChannelPlotOffset() const;
    double             GetArrayPlotOffset() const;
    DisplayTime        GetDisplayTime() const;
    ChannelGeometry    GetChannelGeometry() const;
    double             GetRadius() const;
    double             GetDivergence() const;
    ChannelProfile     GetChannelProfile() const;
    double             GetStandardDeviation() const;
    double             GetSampleDistance() const;
    double             GetSampleVolume() const;
    double             GetSampleArc() const;
    ChannelIntegration GetChannelIntegration() const;
    ToroidalIntegration GetToroidalIntegration() const;
    ToroidalAngleSampling GetToroidalAngleSampling() const;
    double             GetToroidalAngleStart() const;
    double             GetToroidalAngleStop() const;
    double             GetToroidalAngleStride() const;
    TimeSampling       GetTimeSampling() const;
    int                GetTimeStepStart() const;
    int                GetTimeStepStop() const;
    int                GetTimeStepStride() const;
    const doubleVector &GetChannelList() const;
          doubleVector &GetChannelList();
    const doubleVector &GetWallList() const;
          doubleVector &GetWallList();
    int                GetNChannelListArrays() const;
    double             GetChannelListToroidalArrayAngle() const;
    double             GetChannelListToroidalAngle() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string MeshGeometry_ToString(MeshGeometry);
    static bool MeshGeometry_FromString(const std::string &, MeshGeometry &);
protected:
    static std::string MeshGeometry_ToString(int);
public:
    static std::string ArrayConfiguration_ToString(ArrayConfiguration);
    static bool ArrayConfiguration_FromString(const std::string &, ArrayConfiguration &);
protected:
    static std::string ArrayConfiguration_ToString(int);
public:
    static std::string Boundary_ToString(Boundary);
    static bool Boundary_FromString(const std::string &, Boundary &);
protected:
    static std::string Boundary_ToString(int);
public:
    static std::string ChannelProjection_ToString(ChannelProjection);
    static bool ChannelProjection_FromString(const std::string &, ChannelProjection &);
protected:
    static std::string ChannelProjection_ToString(int);
public:
    static std::string ChannelLayoutType_ToString(ChannelLayoutType);
    static bool ChannelLayoutType_FromString(const std::string &, ChannelLayoutType &);
protected:
    static std::string ChannelLayoutType_ToString(int);
public:
    static std::string ArrayAxis_ToString(ArrayAxis);
    static bool ArrayAxis_FromString(const std::string &, ArrayAxis &);
protected:
    static std::string ArrayAxis_ToString(int);
public:
    static std::string ViewGeometry_ToString(ViewGeometry);
    static bool ViewGeometry_FromString(const std::string &, ViewGeometry &);
protected:
    static std::string ViewGeometry_ToString(int);
public:
    static std::string DisplayTime_ToString(DisplayTime);
    static bool DisplayTime_FromString(const std::string &, DisplayTime &);
protected:
    static std::string DisplayTime_ToString(int);
public:
    static std::string ChannelGeometry_ToString(ChannelGeometry);
    static bool ChannelGeometry_FromString(const std::string &, ChannelGeometry &);
protected:
    static std::string ChannelGeometry_ToString(int);
public:
    static std::string ViewDimension_ToString(ViewDimension);
    static bool ViewDimension_FromString(const std::string &, ViewDimension &);
protected:
    static std::string ViewDimension_ToString(int);
public:
    static std::string ChannelProfile_ToString(ChannelProfile);
    static bool ChannelProfile_FromString(const std::string &, ChannelProfile &);
protected:
    static std::string ChannelProfile_ToString(int);
public:
    static std::string ChannelIntegration_ToString(ChannelIntegration);
    static bool ChannelIntegration_FromString(const std::string &, ChannelIntegration &);
protected:
    static std::string ChannelIntegration_ToString(int);
public:
    static std::string ToroidalIntegration_ToString(ToroidalIntegration);
    static bool ToroidalIntegration_FromString(const std::string &, ToroidalIntegration &);
protected:
    static std::string ToroidalIntegration_ToString(int);
public:
    static std::string ToroidalAngleSampling_ToString(ToroidalAngleSampling);
    static bool ToroidalAngleSampling_FromString(const std::string &, ToroidalAngleSampling &);
protected:
    static std::string ToroidalAngleSampling_ToString(int);
public:
    static std::string TimeSampling_ToString(TimeSampling);
    static bool TimeSampling_FromString(const std::string &, TimeSampling &);
protected:
    static std::string TimeSampling_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_meshGeometry = 0,
        ID_arrayConfiguration,
        ID_boundary,
        ID_nArrays,
        ID_toroidalArrayAngle,
        ID_nChannels,
        ID_channelProjection,
        ID_channelLayoutType,
        ID_channelOffset,
        ID_channelAngle,
        ID_nRows,
        ID_rowOffset,
        ID_arrayOrigin,
        ID_arrayAxis,
        ID_poloidalAngleStart,
        ID_poloidalAngleStop,
        ID_poloialAngle,
        ID_poloialRTilt,
        ID_poloialZTilt,
        ID_toroidalAngle,
        ID_flipToroidalAngle,
        ID_viewGeometry,
        ID_viewDimension,
        ID_heightPlotScale,
        ID_channelPlotOffset,
        ID_arrayPlotOffset,
        ID_displayTime,
        ID_channelGeometry,
        ID_radius,
        ID_divergence,
        ID_channelProfile,
        ID_standardDeviation,
        ID_sampleDistance,
        ID_sampleVolume,
        ID_sampleArc,
        ID_channelIntegration,
        ID_toroidalIntegration,
        ID_toroidalAngleSampling,
        ID_toroidalAngleStart,
        ID_toroidalAngleStop,
        ID_toroidalAngleStride,
        ID_timeSampling,
        ID_timeStepStart,
        ID_timeStepStop,
        ID_timeStepStride,
        ID_channelList,
        ID_wallList,
        ID_nChannelListArrays,
        ID_channelListToroidalArrayAngle,
        ID_channelListToroidalAngle,
        ID__LAST
    };

private:
    int          meshGeometry;
    int          arrayConfiguration;
    int          boundary;
    int          nArrays;
    double       toroidalArrayAngle;
    int          nChannels;
    int          channelProjection;
    int          channelLayoutType;
    double       channelOffset;
    double       channelAngle;
    int          nRows;
    double       rowOffset;
    double       arrayOrigin[3];
    int          arrayAxis;
    double       poloidalAngleStart;
    double       poloidalAngleStop;
    double       poloialAngle;
    double       poloialRTilt;
    double       poloialZTilt;
    double       toroidalAngle;
    bool         flipToroidalAngle;
    int          viewGeometry;
    int          viewDimension;
    double       heightPlotScale;
    double       channelPlotOffset;
    double       arrayPlotOffset;
    int          displayTime;
    int          channelGeometry;
    double       radius;
    double       divergence;
    int          channelProfile;
    double       standardDeviation;
    double       sampleDistance;
    double       sampleVolume;
    double       sampleArc;
    int          channelIntegration;
    int          toroidalIntegration;
    int          toroidalAngleSampling;
    double       toroidalAngleStart;
    double       toroidalAngleStop;
    double       toroidalAngleStride;
    int          timeSampling;
    int          timeStepStart;
    int          timeStepStop;
    int          timeStepStride;
    doubleVector channelList;
    doubleVector wallList;
    int          nChannelListArrays;
    double       channelListToroidalArrayAngle;
    double       channelListToroidalAngle;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define LINESAMPLERATTRIBUTES_TMFS "iiiidiiiddidDiddddddbiidddiiddiddddiiidddiiiid*d*idd"

#endif
