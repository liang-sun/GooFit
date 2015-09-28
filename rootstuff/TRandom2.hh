// @(#)root/mathcore:$Id$
// Author: Rene Brun   04/03/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TRandom2
#define ROOT_TRandom2



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRandom2                                                             //
//                                                                      //
// random number generator class (periodicity > 10**26)                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TRandom
#include "TRandom.hh"
#endif

class TRandom2 : public TRandom {

protected:
   unsigned int   fSeed1;  //Random number generator seed 1
   unsigned int   fSeed2;  //Random number generator seed 2

public:
   TRandom2(unsigned int seed=1);
   virtual ~TRandom2();
   virtual  double Rndm(int i=0);
   virtual  void     RndmArray(int n, float *array);
   virtual  void     RndmArray(int n, double *array);
   virtual  void     SetSeed(unsigned int seed=0);

};

#endif
