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

#include <string>
#include "pidx_idx_io.h"
#include "data_handle/PIDX_data_types.h"

typedef std::string String;

int process_count = 1, rank = 0;
//unsigned long long local_box_offset[3];
PIDX_point global_size, local_offset, local_size;

//unsigned long long global_box_size[3] = {0, 0, 0};
//unsigned long long local_box_size[3] = {0, 0, 0};

static void terminate()
{
#if PIDX_HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(-1);
#endif
}

static void terminate_with_error_msg(const char *format, ...)
{
  va_list arg_ptr;
  va_start(arg_ptr, format);
  vfprintf(stderr, format, arg_ptr);
  va_end(arg_ptr);
  terminate();
}

static void init_mpi()
{
#if PIDX_HAVE_MPI
 // if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
 //   terminate_with_error_msg("ERROR: MPI_Init error\n");
  if (MPI_Comm_size(MPI_COMM_WORLD, &process_count) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_size error\n");
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    terminate_with_error_msg("ERROR: MPI_Comm_rank error\n");
#endif
}

VisitIDXIO::DTypes convertType(PIDX_data_type intype){
  
  if(strcmp(intype,INT8) == 0)
    return VisitIDXIO::IDX_INT8;
  else if(strcmp(intype,UINT8) == 0)
    return VisitIDXIO::IDX_UINT8;
  else if(strcmp(intype, INT16)==0)
    return VisitIDXIO::IDX_INT16;
  else if(strcmp(intype, UINT16)==0)
    return VisitIDXIO::IDX_UINT16;
  else if(strcmp(intype, INT32)==0)
    return VisitIDXIO::IDX_INT32;
  else if(strcmp(intype, UINT32)==0)
    return VisitIDXIO::IDX_UINT32;
  else if(strcmp(intype, INT64)==0)
    return VisitIDXIO::IDX_INT64;
  else if(strcmp(intype, UINT64)==0)
    return VisitIDXIO::IDX_UINT64;
  else if(strcmp(intype,FLOAT32)==0)
    return VisitIDXIO::IDX_FLOAT32;
  else if(strcmp(intype, FLOAT64)==0){
    return VisitIDXIO::IDX_FLOAT64;
  }
  
  fprintf(stderr, "Type not found for PIDX\n");
  
}


PIDXIO::~PIDXIO(){
  
}

bool PIDXIO::openDataset(const String filename){
  //init_mpi();
  
  int ret;
  int variable_count;
  PIDX_file file;
  
  //  Creating access
  PIDX_access access;
  PIDX_create_access(&access);
  
  //  PIDX mandatory calls
  ret = PIDX_file_open(filename.c_str(), PIDX_MODE_RDONLY, access, &file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_file_create");
  
  ret = PIDX_get_dims(file, global_size);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_dims");
  
  (global_size[2] > 1) ? dims = 3 : dims = 2;
  
  for(int i=0; i < 3; i++){
    logic_box.p1[i] = 0;
    logic_box.p2[i] = global_size[i]-1;
  }
  
  printf("dims %lld %lld %lld\n", global_size[0],global_size[1],global_size[2]);
  
  ret = PIDX_get_variable_count(file, &variable_count);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_set_variable_count");
  
  int first_tstep = 0, last_tstep = 0;
  ret = PIDX_get_first_tstep(file, &first_tstep);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_first_tstep");
  ret = PIDX_get_last_tstep(file, &last_tstep);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_last_tstep");
  
  for(int i=first_tstep; i <= last_tstep; i++)
    tsteps.push_back(i);
  
  ntimesteps = tsteps.size();
  printf("time size %d\n", tsteps.size());
 /*
  int res[2];
  ret = PIDX_get_resolution(file, res);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_get_first_tstep");
  
  max_resolution = res[1];
  
  printf("max res %d\n", max_resolution);
  */
  printf("%d variables\n", variable_count);
  
  PIDX_variable* variable = (PIDX_variable*)malloc(sizeof(*variable) * variable_count);
  memset(variable, 0, sizeof(*variable) * variable_count);
  
  for (int var = 0; var < variable_count; var++)
  {
    ret = PIDX_get_next_variable(file, &variable[var]);
    if (ret != PIDX_success) terminate_with_error_msg("PIDX_get_next_variable");
    
    int values_per_sample = variable[var]->values_per_sample;
    
    int bits_per_sample = 0;
    ret = PIDX_default_bits_per_datatype(variable[var]->type_name, &bits_per_sample);
    if (ret != PIDX_success)  terminate_with_error_msg("PIDX_default_bytes_per_datatype");
    
    if(rank == 0){
      printf("\t variable %d values per sample %d bits per sample %d \n", var, values_per_sample, bits_per_sample);
    }
    
    VisitIDXIO::Field my_field;
    
    my_field.type = convertType(variable[var]->type_name);
    my_field.isVector = values_per_sample > 1 ? true : false;
    my_field.ncomponents = values_per_sample;
    my_field.name = variable[var]->var_name;
    
    fields.push_back(my_field);
  
  }
  /*
  ret = PIDX_reset_variable_counter(file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_reset_variable_counter");
  
  ret = PIDX_close(file);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close");
  
  ret = PIDX_close_access(access);
  if (ret != PIDX_success)  terminate_with_error_msg("PIDX_close_access");
 */
  return true;

  
}

unsigned char* PIDXIO::getData(const VisitIDXIO::Box box, const int timestate, const char* varname){
  
          printf("GET DATA??\n");
  unsigned char* data = NULL;
  
  return (unsigned char*)data;


}