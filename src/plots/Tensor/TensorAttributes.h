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

#ifndef TENSORATTRIBUTES_H
#define TENSORATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorAttribute.h>

// ****************************************************************************
// Class: TensorAttributes
//
// Purpose:
//    Attributes for the tensor plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class TensorAttributes : public AttributeSubject
{
public:
    TensorAttributes();
    TensorAttributes(const TensorAttributes &obj);
    virtual ~TensorAttributes();

    virtual TensorAttributes& operator = (const TensorAttributes &obj);
    virtual bool operator == (const TensorAttributes &obj) const;
    virtual bool operator != (const TensorAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectTensorColor();
    void SelectColorTableName();

    // Property setting methods
    void SetUseStride(bool useStride_);
    void SetStride(int stride_);
    void SetNTensors(int nTensors_);
    void SetScale(double scale_);
    void SetScaleByMagnitude(bool scaleByMagnitude_);
    void SetAutoScale(bool autoScale_);
    void SetColorByEigenvalues(bool colorByEigenvalues_);
    void SetUseLegend(bool useLegend_);
    void SetTensorColor(const ColorAttribute &tensorColor_);
    void SetColorTableName(const std::string &colorTableName_);

    // Property getting methods
    bool                 GetUseStride() const;
    int                  GetStride() const;
    int                  GetNTensors() const;
    double               GetScale() const;
    bool                 GetScaleByMagnitude() const;
    bool                 GetAutoScale() const;
    bool                 GetColorByEigenvalues() const;
    bool                 GetUseLegend() const;
    const ColorAttribute &GetTensorColor() const;
          ColorAttribute &GetTensorColor();
    const std::string    &GetColorTableName() const;
          std::string    &GetColorTableName();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const TensorAttributes &obj);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_useStride = 0,
        ID_stride,
        ID_nTensors,
        ID_scale,
        ID_scaleByMagnitude,
        ID_autoScale,
        ID_colorByEigenvalues,
        ID_useLegend,
        ID_tensorColor,
        ID_colorTableName
    };

private:
    bool           useStride;
    int            stride;
    int            nTensors;
    double         scale;
    bool           scaleByMagnitude;
    bool           autoScale;
    bool           colorByEigenvalues;
    bool           useLegend;
    ColorAttribute tensorColor;
    std::string    colorTableName;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
