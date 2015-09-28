#include "VariableBinTransform1DPdf.hh"

EXEC_TARGET fptype device_VarBinTransform1D (fptype* evt, fptype* p, unsigned int* indices) {
  // Index structure: nP lim1 lim2 ... 
  int ret = 0;
  //int previousSize = 1; 
  //printf("[%i, %i] Bin Transform: %i %i %f %f\n", THREADIDX, BLOCKIDX, numObservables, previousSize, evt[0], evt[1]); 
  fptype obsValue   = evt[indices[2 + indices[0]]];
  if (obsValue<0) obsValue = -obsValue;
  int numLimits = indices[1];
  for (int i = 0; i < numLimits; ++i) {
    fptype lowerLimit = functorConstants[indices[i+2]];
    if (obsValue < lowerLimit) break;
    ret ++;
  }

  return fptype(ret); 
}

MEM_DEVICE device_function_ptr ptr_to_VarBinTransform1D = device_VarBinTransform1D; 

// Notice that bin sizes and limits can be different, for this purpose, than what's implied by the Variable members. 
__host__ VariableBinTransform1DPdf::VariableBinTransform1DPdf (std::string n, Variable* _x, vector<fptype> binlimits) 
  : GooPdf(_x, n) 
{

  const unsigned int numLimits = binlimits.size(); //Excluding the min & max values for _x
  cIndex = registerConstants(numLimits);
  fptype* host_constants = new fptype[numLimits]; 
  std::vector<unsigned int> pindices;
  pindices.push_back(numLimits); 
  for (unsigned int i = 0; i < numLimits; ++i) {
    pindices.push_back(cIndex + i); 

    host_constants[i] = binlimits[i]; // cIndex will be accounted for by offset in memcpy
  }

  MEMCPY_TO_SYMBOL(functorConstants, host_constants, numLimits*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice); 
  delete[] host_constants; 

  GET_FUNCTION_ADDR(ptr_to_VarBinTransform1D);
  initialise(pindices); 
}


