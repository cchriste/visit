/***************************************************
 ** ViSUS Visualization Project                    **
 ** Copyright (c) 2010 University of Utah          **
 ** Scientific Computing and Imaging Institute     **
 ** 72 S Central Campus Drive, Room 3750           **
 ** Salt Lake City, UT 84112                       **
 **                                                **
 ** For information about this project see:        **
 ** http://www.pascucci.org/visus/                 **
 **                                                **
 **      or contact: pascucci@sci.utah.edu         **
 **                                                **
 ****************************************************/

#ifndef _visus_simpleio_h
#define _visus_simpleio_h

#include <string>
#include <vector>
#include <cassert>
#include "visus_simpleio_types.h"

using namespace VisusSimpleIO;

class DatasetImpl;
class AccessImpl;

// TODO generalize end extend
// Query at full resolution only

class SimpleIO{
    
public:
    
    SimpleIO(){};
    
    bool openDataset(const std::string filename);
    
    inline int getDimension(){
        return dims;
    }
    
    inline int getNTimesteps(){
        return ntimesteps;
    }
    
    inline int getMaxResolution(){
        return max_resolution;
    }
    
    unsigned char* getData(const SimpleBox box, const int timestate, const char* varname, bool reverseEndianess);
    
    inline std::vector<double> getTimes(){
        return tsteps;
    }
    
    inline std::vector<SimpleField> getFields(){
        return fields;
    }
    
    inline SimpleBox getLogicBox(){
        return logic_box;
    }
    
    inline const double* getLogicToPhysic(){
        return &logic_to_physic[0];
    }
    
    inline SimpleField getCurrField(){
        return curr_field;
    }
  
    inline bool isCompressed(){
      return compressed_dataset;
    }
    
    ~SimpleIO();
    
private:
    int dims;
    int ntimesteps;
    int max_resolution;
    std::vector<double> tsteps;
    std::string dataset_url;
    std::vector<SimpleField> fields;
    SimpleBox logic_box;
    SimpleField curr_field;
    double logic_to_physic[16];
    DatasetImpl* datasetImpl;
    bool compressed_dataset;
};


#endif
