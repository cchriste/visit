/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
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
// Creation:   Thu Mar 8 15:21:18 PST 2007
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
    enum TwoDAmount
    {
        Area,
        RevolvedVolume
    };
    enum BasedOn
    {
        ManyVarsForSingleZone,
        ManyZonesForSingleVar
    };
    enum BinContribution
    {
        Auto,
        Frequency,
        Weighted
    };

    HistogramAttributes();
    HistogramAttributes(const HistogramAttributes &obj);
    virtual ~HistogramAttributes();

    virtual HistogramAttributes& operator = (const HistogramAttributes &obj);
    virtual bool operator == (const HistogramAttributes &obj) const;
    virtual bool operator != (const HistogramAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectColor();

    // Property setting methods
    void SetBasedOn(BasedOn basedOn_);
    void SetHistogramType(BinContribution histogramType_);
    void SetTwoDAmount(TwoDAmount twoDAmount_);
    void SetSpecifyRange(bool specifyRange_);
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

    // Property getting methods
    BasedOn              GetBasedOn() const;
    BinContribution      GetHistogramType() const;
    TwoDAmount           GetTwoDAmount() const;
    bool                 GetSpecifyRange() const;
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

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string OutputType_ToString(OutputType);
    static bool OutputType_FromString(const std::string &, OutputType &);
protected:
    static std::string OutputType_ToString(int);
public:
    static std::string TwoDAmount_ToString(TwoDAmount);
    static bool TwoDAmount_FromString(const std::string &, TwoDAmount &);
protected:
    static std::string TwoDAmount_ToString(int);
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

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const HistogramAttributes &) const;
    virtual bool VarChangeRequiresReset(void);
private:
    int            basedOn;
    int            histogramType;
    int            twoDAmount;
    bool           specifyRange;
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
};

#endif
