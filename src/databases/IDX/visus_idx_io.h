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

#ifndef _visus_idx_io_h
#define _visus_idx_io_h

#include <string>
#include <vector>
#include <cassert>
#include "visit_idx_io.h"
#include "visit_idx_io_types.h"

using namespace VisitIDXIO;

class DatasetImpl;
class AccessImpl;

// TODO generalize end extend
// Query at full resolution only

class VisusIDXIO : public IDX_IO{
    
public:
    
    VisusIDXIO(){};
    
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
    
    unsigned char* getData(const Box box, const int timestate, const char* varname);
    
    inline std::vector<double> getTimes(){
        return tsteps;
    }
    
    inline std::vector<Field> getFields(){
        return fields;
    }
    
    inline Box getLogicBox(){
        return logic_box;
    }
    
    inline const double* getLogicToPhysic(){
        return &logic_to_physic[0];
    }
    
    inline Field getCurrField(){
        return curr_field;
    }
  
    inline bool isCompressed(){
      return compressed_dataset;
    }
    
    virtual ~VisusIDXIO();
  
private:
  DatasetImpl* datasetImpl;
  
};


#endif
