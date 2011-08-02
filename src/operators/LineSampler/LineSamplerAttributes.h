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
    enum CoordinateSystem
    {
        Cartesian,
        Cylindrical
    };
    enum ArrayProjection
    {
        Parallel,
        Divergent,
        Grid
    };
    enum ArrayAxis
    {
        R,
        Z
    };
    enum BeamShape
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
    enum BeamType
    {
        TopHat,
        Gaussian
    };
    enum ChannelSampling
    {
        SingleChannelSampling,
        SummationChannelSampling
    };
    enum ToroidalSampling
    {
        SingleToroidalSampling,
        SummationToroidalSampling
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

    // Property setting methods
    void SetCoordinateSystem(CoordinateSystem coordinateSystem_);
    void SetNArrays(int nArrays_);
    void SetNChannels(int nChannels_);
    void SetToroialArrayAngle(double toroialArrayAngle_);
    void SetArrayProjection(ArrayProjection arrayProjection_);
    void SetChannelOffset(double channelOffset_);
    void SetChannelAngle(double channelAngle_);
    void SetNRows(int nRows_);
    void SetRowOffset(double rowOffset_);
    void SetArrayOrigin(const double *arrayOrigin_);
    void SetArrayAxis(ArrayAxis arrayAxis_);
    void SetPoloialAngle(double poloialAngle_);
    void SetPoloialRTilt(double poloialRTilt_);
    void SetPoloialZTilt(double poloialZTilt_);
    void SetToroialAngle(double toroialAngle_);
    void SetViewDimension(ViewDimension viewDimension_);
    void SetHeightPlotScale(double heightPlotScale_);
    void SetChannelPlotOffset(double channelPlotOffset_);
    void SetArrayPlotOffset(double arrayPlotOffset_);
    void SetTimePlotScale(double timePlotScale_);
    void SetBeamShape(BeamShape beamShape_);
    void SetRadius(double radius_);
    void SetDivergence(double divergence_);
    void SetBeamType(BeamType beamType_);
    void SetStandardDeviation(double standardDeviation_);
    void SetSampleDistance(double sampleDistance_);
    void SetSampleArc(double sampleArc_);
    void SetChannelSampling(ChannelSampling channelSampling_);
    void SetToroidalSampling(ToroidalSampling toroidalSampling_);
    void SetToroidalSamplingAngle(double toroidalSamplingAngle_);
    void SetTimeSampling(TimeSampling timeSampling_);
    void SetTimeStepStart(int timeStepStart_);
    void SetTimeStepStop(int timeStepStop_);
    void SetTimeStepStride(int timeStepStride_);

    // Property getting methods
    CoordinateSystem GetCoordinateSystem() const;
    int          GetNArrays() const;
    int          GetNChannels() const;
    double       GetToroialArrayAngle() const;
    ArrayProjection GetArrayProjection() const;
    double       GetChannelOffset() const;
    double       GetChannelAngle() const;
    int          GetNRows() const;
    double       GetRowOffset() const;
    const double *GetArrayOrigin() const;
          double *GetArrayOrigin();
    ArrayAxis    GetArrayAxis() const;
    double       GetPoloialAngle() const;
    double       GetPoloialRTilt() const;
    double       GetPoloialZTilt() const;
    double       GetToroialAngle() const;
    ViewDimension GetViewDimension() const;
    double       GetHeightPlotScale() const;
    double       GetChannelPlotOffset() const;
    double       GetArrayPlotOffset() const;
    double       GetTimePlotScale() const;
    BeamShape    GetBeamShape() const;
    double       GetRadius() const;
    double       GetDivergence() const;
    BeamType     GetBeamType() const;
    double       GetStandardDeviation() const;
    double       GetSampleDistance() const;
    double       GetSampleArc() const;
    ChannelSampling GetChannelSampling() const;
    ToroidalSampling GetToroidalSampling() const;
    double       GetToroidalSamplingAngle() const;
    TimeSampling GetTimeSampling() const;
    int          GetTimeStepStart() const;
    int          GetTimeStepStop() const;
    int          GetTimeStepStride() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string CoordinateSystem_ToString(CoordinateSystem);
    static bool CoordinateSystem_FromString(const std::string &, CoordinateSystem &);
protected:
    static std::string CoordinateSystem_ToString(int);
public:
    static std::string ArrayProjection_ToString(ArrayProjection);
    static bool ArrayProjection_FromString(const std::string &, ArrayProjection &);
protected:
    static std::string ArrayProjection_ToString(int);
public:
    static std::string ArrayAxis_ToString(ArrayAxis);
    static bool ArrayAxis_FromString(const std::string &, ArrayAxis &);
protected:
    static std::string ArrayAxis_ToString(int);
public:
    static std::string BeamShape_ToString(BeamShape);
    static bool BeamShape_FromString(const std::string &, BeamShape &);
protected:
    static std::string BeamShape_ToString(int);
public:
    static std::string ViewDimension_ToString(ViewDimension);
    static bool ViewDimension_FromString(const std::string &, ViewDimension &);
protected:
    static std::string ViewDimension_ToString(int);
public:
    static std::string BeamType_ToString(BeamType);
    static bool BeamType_FromString(const std::string &, BeamType &);
protected:
    static std::string BeamType_ToString(int);
public:
    static std::string ChannelSampling_ToString(ChannelSampling);
    static bool ChannelSampling_FromString(const std::string &, ChannelSampling &);
protected:
    static std::string ChannelSampling_ToString(int);
public:
    static std::string ToroidalSampling_ToString(ToroidalSampling);
    static bool ToroidalSampling_FromString(const std::string &, ToroidalSampling &);
protected:
    static std::string ToroidalSampling_ToString(int);
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
        ID_coordinateSystem = 0,
        ID_nArrays,
        ID_nChannels,
        ID_toroialArrayAngle,
        ID_arrayProjection,
        ID_channelOffset,
        ID_channelAngle,
        ID_nRows,
        ID_rowOffset,
        ID_arrayOrigin,
        ID_arrayAxis,
        ID_poloialAngle,
        ID_poloialRTilt,
        ID_poloialZTilt,
        ID_toroialAngle,
        ID_viewDimension,
        ID_heightPlotScale,
        ID_channelPlotOffset,
        ID_arrayPlotOffset,
        ID_timePlotScale,
        ID_beamShape,
        ID_radius,
        ID_divergence,
        ID_beamType,
        ID_standardDeviation,
        ID_sampleDistance,
        ID_sampleArc,
        ID_channelSampling,
        ID_toroidalSampling,
        ID_toroidalSamplingAngle,
        ID_timeSampling,
        ID_timeStepStart,
        ID_timeStepStop,
        ID_timeStepStride,
        ID__LAST
    };

private:
    int    coordinateSystem;
    int    nArrays;
    int    nChannels;
    double toroialArrayAngle;
    int    arrayProjection;
    double channelOffset;
    double channelAngle;
    int    nRows;
    double rowOffset;
    double arrayOrigin[3];
    int    arrayAxis;
    double poloialAngle;
    double poloialRTilt;
    double poloialZTilt;
    double toroialAngle;
    int    viewDimension;
    double heightPlotScale;
    double channelPlotOffset;
    double arrayPlotOffset;
    double timePlotScale;
    int    beamShape;
    double radius;
    double divergence;
    int    beamType;
    double standardDeviation;
    double sampleDistance;
    double sampleArc;
    int    channelSampling;
    int    toroidalSampling;
    double toroidalSamplingAngle;
    int    timeSampling;
    int    timeStepStart;
    int    timeStepStop;
    int    timeStepStride;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define LINESAMPLERATTRIBUTES_TMFS "iiididdidDiddddiddddiddidddiidiiii"

#endif
