
/* 
  This file contains utility functions to handle Uintah metadata
*/

#ifndef UINTAH_IDX_UTILS_H
#define UINTAH_IDX_UTILS_H

#include "LevelInfo.h"
#include <vtkSmartPointer.h>
#include <vtkXMLDataParser.h>
#include <sstream>

extern void ups_parse_vector(vtkXMLDataElement *el, int* vec, int dim);
extern void ups_parse_vector(vtkXMLDataElement *el, double* vec, int dim);

extern void parse_ups(vtkSmartPointer<vtkXMLDataParser> parser, LevelInfo& levelInfo, int dim, bool use_extracells);

static inline int        cint   (std::string s) {int    value;std::istringstream iss(s);iss>>value;return value;}
static inline float      cfloat (std::string s) {float  value;std::istringstream iss(s);iss>>value;return value;}
static inline double     cdouble(std::string s) {double value;std::istringstream iss(s);iss>>value;return value;}
// trim from start
static inline std::string &ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
  return ltrim(rtrim(s));
}

#endif //UINTAH_IDX_UTILS_H
