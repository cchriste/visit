/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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

#ifndef HISTOGRAMATTRIBUTES_H
#define HISTOGRAMATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

#include <ColorAttribute.h>
#include <PickAttributes.h>

// ****************************************************************************
// Class: HistogramAttributes
//
// Purpose:
//    Attributes for Histogram Plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class HistogramAttributes : public AttributeSubject
{
public:
    enum OutputType
    {
        Curve,
        Block
    };
    enum BasedOn
    {
        ManyVarsForSingleZone,
        ManyZonesForSingleVar
    };
    enum BinContribution
    {
        Frequency,
        Weighted,
        Variable
    };
    enum LimitsMode
    {
        OriginalData,
        CurrentPlot
    };
    enum DataScale
    {
        Linear,
        Log,
        SquareRoot
    };

    // These constructors are for objects of this class
    HistogramAttributes();
    HistogramAttributes(const HistogramAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    HistogramAttributes(private_tmfs_t tmfs);
    HistogramAttributes(const HistogramAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~HistogramAttributes();

    virtual HistogramAttributes& operator = (const HistogramAttributes &obj);
    virtual bool operator == (const HistogramAttributes &obj) const;
    virtual bool operator != (const HistogramAttributes &obj) const;
private:
    void Init();
    void Copy(const HistogramAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectWeightVariable();
    void SelectColor();

    // Property setting methods
    void SetBasedOn(BasedOn basedOn_);
    void SetHistogramType(BinContribution histogramType_);
    void SetWeightVariable(const std::string &weightVariable_);
    void SetLimitsMode(LimitsMode limitsMode_);
    void SetMinFlag(bool minFlag_);
    void SetMaxFlag(bool maxFlag_);
    void SetMin(double min_);
    void SetMax(double max_);
    void SetNumBins(int numBins_);
    void SetDomain(int domain_);
    void SetZone(int zone_);
    void SetUseBinWidths(bool useBinWidths_);
    void SetOutputType(OutputType outputType_);
    void SetLineStyle(int lineStyle_);
    void SetLineWidth(int lineWidth_);
    void SetColor(const ColorAttribute &color_);
    void SetDataScale(DataScale dataScale_);
    void SetBinScale(DataScale binScale_);
    void SetNormalizeHistogram(bool normalizeHistogram_);
    void SetComputeAsCDF(bool computeAsCDF_);

    // Property getting methods
    BasedOn              GetBasedOn() const;
    BinContribution      GetHistogramType() const;
    const std::string    &GetWeightVariable() const;
          std::string    &GetWeightVariable();
    LimitsMode           GetLimitsMode() const;
    bool                 GetMinFlag() const;
    bool                 GetMaxFlag() const;
    double               GetMin() const;
    double               GetMax() const;
    int                  GetNumBins() const;
    int                  GetDomain() const;
    int                  GetZone() const;
    bool                 GetUseBinWidths() const;
    OutputType           GetOutputType() const;
    int                  GetLineStyle() const;
    int                  GetLineWidth() const;
    const ColorAttribute &GetColor() const;
          ColorAttribute &GetColor();
    DataScale            GetDataScale() const;
    DataScale            GetBinScale() const;
    bool                 GetNormalizeHistogram() const;
    bool                 GetComputeAsCDF() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string OutputType_ToString(OutputType);
    static bool OutputType_FromString(const std::string &, OutputType &);
protected:
    static std::string OutputType_ToString(int);
public:
    static std::string BasedOn_ToString(BasedOn);
    static bool BasedOn_FromString(const std::string &, BasedOn &);
protected:
    static std::string BasedOn_ToString(int);
public:
    static std::string BinContribution_ToString(BinContribution);
    static bool BinContribution_FromString(const std::string &, BinContribution &);
protected:
    static std::string BinContribution_ToString(int);
public:
    static std::string LimitsMode_ToString(LimitsMode);
    static bool LimitsMode_FromString(const std::string &, LimitsMode &);
protected:
    static std::string LimitsMode_ToString(int);
public:
    static std::string DataScale_ToString(DataScale);
    static bool DataScale_FromString(const std::string &, DataScale &);
protected:
    static std::string DataScale_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const HistogramAttributes &) const;
    virtual bool VarChangeRequiresReset(void);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_basedOn = 0,
        ID_histogramType,
        ID_weightVariable,
        ID_limitsMode,
        ID_minFlag,
        ID_maxFlag,
        ID_min,
        ID_max,
        ID_numBins,
        ID_domain,
        ID_zone,
        ID_useBinWidths,
        ID_outputType,
        ID_lineStyle,
        ID_lineWidth,
        ID_color,
        ID_dataScale,
        ID_binScale,
        ID_normalizeHistogram,
        ID_computeAsCDF,
        ID__LAST
    };

private:
    int            basedOn;
    int            histogramType;
    std::string    weightVariable;
    int            limitsMode;
    bool           minFlag;
    bool           maxFlag;
    double         min;
    double         max;
    int            numBins;
    int            domain;
    int            zone;
    bool           useBinWidths;
    int            outputType;
    int            lineStyle;
    int            lineWidth;
    ColorAttribute color;
    int            dataScale;
    int            binScale;
    bool           normalizeHistogram;
    bool           computeAsCDF;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define HISTOGRAMATTRIBUTES_TMFS "iisibbddiiibiiiaiibb"

#endif
