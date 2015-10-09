/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2008 Scientific Computing and Imaging Institute,
   University of Utah.

   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/

/**
  \file    TransferFunction2D.cpp
  \author  Jens Krueger
           SCI Institute
           University of Utah
  \version 1.0
  \date    September 2008
*/

#include <memory.h>
#include "TransferFunction2D.h"

#include "Controller/Controller.h"

using namespace std;

TransferFunction2D::TransferFunction2D() :
  m_iSize(0,0),
  m_pColorData(NULL),
  m_pPixelData(NULL),
  m_pRCanvas(NULL),
  m_vValueBBox(0,0),
  m_bUseCachedData(false)
{
  Resize(m_iSize);
}

TransferFunction2D::TransferFunction2D(const std::string& filename):
  m_iSize(0,0),
  m_pColorData(NULL),
  m_pPixelData(NULL),
  m_pRCanvas(NULL),
  m_bUseCachedData(false)
{
  Load(filename);
}

TransferFunction2D::TransferFunction2D(const VECTOR2<size_t>& iSize):
  m_iSize(iSize),
  m_pColorData(NULL),
  m_pPixelData(NULL),
  m_pRCanvas(NULL),
  m_bUseCachedData(false)
{
  Resize(m_iSize);
}


void TransferFunction2D::DeleteCanvasData()
{
  delete m_pColorData;
  delete [] m_pPixelData;
  delete m_pRCanvas;
}

TransferFunction2D::~TransferFunction2D(void)
{
  DeleteCanvasData();
}

void TransferFunction2D::Resize(const VECTOR2<size_t>& iSize) {
  m_iSize = iSize;
  m_Trans1D.Resize(iSize.x);
  m_Trans1D.Clear();

  DeleteCanvasData();
}

void TransferFunction2D::Resample(const VECTOR2<size_t>& iSize) {
  m_iSize = iSize;
  m_Trans1D.Resample(iSize.x);
}

bool TransferFunction2D::Load(const std::string& filename, const VECTOR2<size_t>& vTargetSize) {
  ifstream file(filename.c_str());

  if (!file.is_open()) return false;

  m_iSize = vTargetSize;

  // ignore the size in the file (read it but never use it again)
  VECTOR2<size_t> vSizeInFile;
  file >> vSizeInFile.x >> vSizeInFile.y;

  // load 1D Trans
  m_Trans1D.Load(file, vTargetSize.x);

  // load swatch count
  UINT32 iSwatchCount;
  file >> iSwatchCount;
  m_Swatches.resize(iSwatchCount);

  // load Swatches
  for (size_t i = 0;i<m_Swatches.size();i++) m_Swatches[i].Load(file);

  file.close();

  return true;
}


bool TransferFunction2D::Load(const std::string& filename) {
  ifstream file(filename.c_str());

  if (!file.is_open()) return false;

  // load size
  file >> m_iSize.x >> m_iSize.y;

  // load 1D Trans
  m_Trans1D.Load(file, m_iSize.x);

  // load swatch count
  UINT32 iSwatchCount;
  file >> iSwatchCount;
  m_Swatches.resize(iSwatchCount);

  // load Swatches
  for (size_t i = 0;i<m_Swatches.size();i++) m_Swatches[i].Load(file);

  file.close();

  return true;
}

bool TransferFunction2D::Save(const std::string& filename) const {
  ofstream file(filename.c_str());

  if (!file.is_open()) return false;

  // save size
  file << m_iSize.x << " " << m_iSize.y << endl;

  // save 1D Trans
  m_Trans1D.Save(file);

  // save swatch count
  file << m_Swatches.size() << endl;

  // save Swatches
  for (size_t i = 0;i<m_Swatches.size();i++) m_Swatches[i].Save(file);

  file.close();

  return true;
}

void TransferFunction2D::GetByteArray(unsigned char** pcData) {
  if (*pcData == NULL) *pcData = new unsigned char[m_iSize.area()*4];

  size_t iSize = m_iSize.area();
  unsigned char *pcSourceDataIterator = RenderTransferFunction8Bit();
  unsigned char *pcDataIterator = *pcData;
  memcpy(pcDataIterator, pcSourceDataIterator, iSize*4);
  for (size_t i = 0;i<iSize;i++) {
    unsigned char r = *(pcDataIterator+2);
    unsigned char b = *(pcDataIterator+0);

    *(pcDataIterator+0) = r;
    *(pcDataIterator+2) = b;

    pcDataIterator+=4;
  }
}


void TransferFunction2D::GetByteArray(unsigned char** pcData,
                                      unsigned char cUsedRange) {
  if (*pcData == NULL) *pcData = new unsigned char[m_iSize.area()*4];

  float fScale = 255.0f/float(cUsedRange);

  size_t iSize = m_iSize.area();
  unsigned char *pcSourceDataIterator = RenderTransferFunction8Bit();
  unsigned char *pcDataIterator = *pcData;
  memcpy(pcDataIterator, pcSourceDataIterator, iSize*4);
  for (size_t i = 0;i<iSize;i++) {
    unsigned char r = *(pcDataIterator+2);
    unsigned char g = *(pcDataIterator+1);
    unsigned char b = *(pcDataIterator+0);
    unsigned char a = *(pcDataIterator+3);

    *(pcDataIterator+0) = (unsigned char)(float(r)*fScale);
    *(pcDataIterator+1) = (unsigned char)(float(g)*fScale);
    *(pcDataIterator+2) = (unsigned char)(float(b)*fScale);
    *(pcDataIterator+3) = (unsigned char)(float(a)*fScale);

    pcDataIterator+=4;
  }
}

void TransferFunction2D::GetShortArray(unsigned short** psData,
                                       unsigned short sUsedRange) {
  if (*psData == NULL) *psData = new unsigned short[m_iSize.area()*4];

  RenderTransferFunction();
  unsigned short *psDataIterator = *psData;
  FLOATVECTOR4  *piSourceIterator = m_pColorData->GetDataPointer();
  for (size_t i = 0;i<m_pColorData->GetSize().area();i++) {
    *psDataIterator++ = (unsigned short)((*piSourceIterator)[0]*sUsedRange);
    *psDataIterator++ = (unsigned short)((*piSourceIterator)[1]*sUsedRange);
    *psDataIterator++ = (unsigned short)((*piSourceIterator)[2]*sUsedRange);
    *psDataIterator++ = (unsigned short)((*piSourceIterator)[3]*sUsedRange);
    piSourceIterator++;
  }
}

void TransferFunction2D::GetFloatArray(float** pfData) {
  if (*pfData == NULL) *pfData = new float[4*m_iSize.area()];

  RenderTransferFunction();
  memcpy(*pfData, m_pColorData->GetDataPointer(), 4*sizeof(float)*m_iSize.area());
}

INTVECTOR2 TransferFunction2D::Normalized2Offscreen(FLOATVECTOR2 vfCoord, VECTOR2<size_t> iSize) const {
  return INTVECTOR2(int(vfCoord.x*int(iSize.x)),
                    int(vfCoord.y*int(iSize.y)));
}

unsigned char* TransferFunction2D::RenderTransferFunction8Bit() {
#ifndef TUVOK_NO_QT
  VECTOR2<size_t> vRS = GetRenderSize();
  if (m_pColorData == NULL ) m_pColorData = new ColorData2D(m_iSize);
  if (m_pPixelData == NULL ) m_pPixelData = new unsigned char[4*m_iSize.area()];
  if (m_pRCanvas == NULL )   m_pRCanvas   = new QImage(int(vRS.x), int(vRS.y), QImage::Format_ARGB32);

  if (m_pPixelData != NULL && m_bUseCachedData) return m_pPixelData;

  m_pRCanvas->fill(0);

  // render 1D trans 
  QRect imageRect(0, 0, int(vRS.x), int(vRS.y));
  m_Painter.begin(m_pRCanvas);
  m_Painter.drawImage(imageRect,m_Trans1DImage);

  // render swatches
  QPen noBorderPen(Qt::NoPen);
  m_Painter.setPen(noBorderPen);
  for (size_t i = 0;i<m_Swatches.size();i++) {
    TFPolygon& currentSwatch = m_Swatches[i];
    
    std::vector<QPoint> pointList(currentSwatch.pPoints.size());
    for (size_t j = 0;j<currentSwatch.pPoints.size();j++) {    
      INTVECTOR2 vPixelPos = Normalized2Offscreen(currentSwatch.pPoints[j],vRS);
      pointList[j] = QPoint(vPixelPos.x, vPixelPos.y);
    }

    INTVECTOR2 vPixelPos0 = Normalized2Offscreen(currentSwatch.pGradientCoords[0],vRS), 
               vPixelPos1 = Normalized2Offscreen(currentSwatch.pGradientCoords[1],vRS); 

    QGradient* pGradientBrush;
    if (currentSwatch.bRadial) {
      double r = sqrt( pow(double(vPixelPos0.x-vPixelPos1.x),2.0) + pow(double(vPixelPos0.y-vPixelPos1.y),2.0));
      pGradientBrush = new QRadialGradient(vPixelPos0.x, vPixelPos0.y, r);
    } else {
      pGradientBrush = new QLinearGradient(vPixelPos0.x, vPixelPos0.y, vPixelPos1.x, vPixelPos1.y);
    }
    
    for (size_t j = 0;j<currentSwatch.pGradientStops.size();j++) {      
      pGradientBrush->setColorAt(currentSwatch.pGradientStops[j].first, 
                   QColor(int(currentSwatch.pGradientStops[j].second[0]*255),
                          int(currentSwatch.pGradientStops[j].second[1]*255),
                          int(currentSwatch.pGradientStops[j].second[2]*255),
                          int(currentSwatch.pGradientStops[j].second[3]*255)));
    }

    m_Painter.setBrush(*pGradientBrush);
    if(!pointList.empty()) {
      m_Painter.drawPolygon(&pointList[0], int(currentSwatch.pPoints.size()));
    }
    delete pGradientBrush;
  }
  m_Painter.end();


  memcpy(m_pPixelData, m_pRCanvas->scaled(int(m_iSize.x), int(m_iSize.y)).bits(), 4*m_iSize.area());
  m_bUseCachedData = true;

  return m_pPixelData;
#else
  T_ERROR("Cannot render transfer functions without Qt.");
  return NULL;
#endif
}

ColorData2D* TransferFunction2D::RenderTransferFunction() {
  
  unsigned char* pPixelData = RenderTransferFunction8Bit();

  FLOATVECTOR4* p = (FLOATVECTOR4*)(m_pColorData->GetDataPointer());
  for (size_t i = 0;i<m_pColorData->GetSize().area();i++) {
    p[i] = FLOATVECTOR4(pPixelData[4*i+2]/255.0f,
                        pPixelData[4*i+1]/255.0f,
                        pPixelData[4*i+0]/255.0f,
                        pPixelData[4*i+3]/255.0f);
  }

  return m_pColorData;
}

void TransferFunction2D::ComputeNonZeroLimits() {   
  unsigned char* pPixelData = RenderTransferFunction8Bit();

  m_vValueBBox    = UINT64VECTOR4(UINT64(m_iSize.x),0,
                                  UINT64(m_iSize.y),0);

  size_t i = 3;
  for (size_t y = 0;y<m_iSize.y;y++) {
    for (size_t x = 0;x<m_iSize.x;x++) {
      if (pPixelData[i] != 0) {
        m_vValueBBox.x = MIN(m_vValueBBox.x, x);
        m_vValueBBox.y = MAX(m_vValueBBox.y, x);

        m_vValueBBox.z = MIN(m_vValueBBox.z, y);
        m_vValueBBox.w = MAX(m_vValueBBox.w, y);
      }
      i+=4;
    }
  }
}

void TransferFunction2D::Update1DTrans(const TransferFunction1D* p1DTrans) {
#ifndef TUVOK_NO_QT
  m_Trans1D = TransferFunction1D(*p1DTrans);

  size_t iSize = min<size_t>(m_iSize.x,  m_Trans1D.GetSize());

  m_Trans1DImage = QImage(int(iSize), 1, QImage::Format_ARGB32);
  for (size_t i = 0;i<iSize;i++) {
    m_Trans1DImage.setPixel(int(i),0,qRgba(int(m_Trans1D.vColorData[i][0]*255),
                                           int(m_Trans1D.vColorData[i][1]*255),
                                           int(m_Trans1D.vColorData[i][2]*255),
                                           int(m_Trans1D.vColorData[i][3]*255)));
  }
#else
  T_ERROR("Unsupported without Qt.");
#endif
}

// ***************************************************************************

void TFPolygon::Load(ifstream& file) {
  UINT32 iSize;
  file >> bRadial;
  file >> iSize;
  pPoints.resize(iSize);

  for(size_t i=0;i<pPoints.size();++i){
    for(size_t j=0;j<2;++j){
      file >> pPoints[i][j];
    }
  }

  file >> pGradientCoords[0][0] >> pGradientCoords[0][1];
  file >> pGradientCoords[1][0] >> pGradientCoords[1][1];

  file >> iSize;
  pGradientStops.resize(iSize);
  for(size_t i=0;i<pGradientStops.size();++i){
    file >> pGradientStops[i].first;
    for(size_t j=0;j<4;++j){
      file >> pGradientStops[i].second[j];
    }
  }

}

void TFPolygon::Save(ofstream& file) const {
  file << bRadial << endl;
  file << UINT32(pPoints.size()) << endl;

  for(size_t i=0;i<pPoints.size();++i){
    for(size_t j=0;j<2;++j){
      file << pPoints[i][j] << " ";
    }
    file << endl;
  }

  file << pGradientCoords[0][0] << " " << pGradientCoords[0][1] << " ";
  file << pGradientCoords[1][0] << " " << pGradientCoords[1][1];
  file << endl;
  file << pGradientStops.size() << endl;

  for(size_t i=0;i<pGradientStops.size();++i){
    file << pGradientStops[i].first << "  ";
    for(size_t j=0;j<4;++j){
      file << pGradientStops[i].second[j] << " ";
    }
    file << endl;
  }  
}
