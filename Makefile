#------------------------------------------------------------------------------
CXX=nvcc
LD=g++ 
OutPutOpt = -o

#CXXFLAGS     = -dc -std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -D__STRICT_ANSI__ -O3 #-g -G 
CXXFLAGS     = -dc -std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -O3 #-g -G 
DEFINEFLAGS=-DDUMMY=dummy 
CUDALIBDIR=lib64

UNAME=$(shell uname)
ifeq ($(UNAME), Darwin)
CUDALIBDIR=lib
CXXFLAGS+=-m64
endif

EIGENINCDIR=$(HOME)/MinuitThrust/eigen3

ifneq ($(CUDAPRINT),)
DEFINEFLAGS += -DCUDAPRINT=yes
endif 

ifneq ($(PRINTCALLS),)
DEFINEFLAGS += -DPRINTCALLS=yes
endif 

ifneq ($(PROFILE),)
DEFINEFLAGS += -DPROFILING=yes
endif 

ifeq ($(TARGET_OMP),)
CXXFLAGS += -arch=sm_20
else
DEFINEFLAGS += -fno-inline -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_BACKEND_OMP
LIBS += -lgomp
endif 

ifeq ($(CUDALOCATION), )
CUDALOCATION = /usr/local/cuda/
endif
CUDAHEADERS = $(CUDALOCATION)/include/

PWD = $(shell /bin/pwd)
SRCDIR = $(PWD)/PDFs

INCLUDES += -I$(SRCDIR) -I$(PWD) -I$(CUDAHEADERS) -I$(PWD)/rootstuff -I$(PWD)/MCBooster
INCLUDES += -I$(EIGENINCDIR)
LIBS += -L$(CUDALOCATION)/$(CUDALIBDIR) -lcudart -L$(PWD)/rootstuff -lRootUtils 

# GooPdf must be first in CUDAglob, as it defines global variables.
FUNCTORLIST    = $(SRCDIR)/GooPdf.cu 
FUNCTORLIST    += $(filter-out $(SRCDIR)/GooPdf.cu, $(wildcard $(SRCDIR)/*Pdf.cu))
FUNCTORLIST   += $(wildcard $(SRCDIR)/*Aux.cu)
HEADERLIST     = $(patsubst %.cu,%.hh,$(SRCFILES))
WRKFUNCTORLIST = $(patsubst $(SRCDIR)/%.cu,wrkdir/%.cu,$(FUNCTORLIST))
#NB, the above are used in the SRCDIR Makefile.

THRUSTO		= wrkdir/Variable.o wrkdir/FitManager.o wrkdir/GooPdfCUDA.o wrkdir/Faddeeva.o wrkdir/FitControl.o wrkdir/PdfBase.o wrkdir/DataSet.o wrkdir/BinnedDataSet.o wrkdir/UnbinnedDataSet.o wrkdir/FunctorWriter.o 
LINKOBJ         = link.o
ROOTRIPDIR	= $(PWD)/rootstuff
ROOTRIPOBJS	= $(ROOTRIPDIR)/TMinuit.o $(ROOTRIPDIR)/TRandom.o $(ROOTRIPDIR)/TRandom3.o $(ROOTRIPDIR)/TRandom2.o
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 
GOOFITLIB	= libGooFit.a

.SUFFIXES: 
.PHONY:		goofit clean 

goofit:		$(THRUSTO) $(ROOTUTILLIB) 
		@echo "Built GooFit objects" 
		$(CXX) -dlink -arch=sm_20 -o $(LINKOBJ) $(THRUSTO)
		$(CXX) -lib $(LINKOBJ) $(THRUSTO) -o $(GOOFITLIB)

# One rule for GooFit objects.
wrkdir/%.o:	%.cc %.hh 
		@mkdir -p wrkdir 
		$(CXX) $(INCLUDES) $(CXXFLAGS) $(DEFINEFLAGS)  -o $@ -c $<

# A different rule for user-level objects. Notice ROOT_INCLUDES. 
%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS)  -o $@ -c $<

# Still a third rule for the ROOT objects - these have their own Makefile. 
$(ROOTRIPDIR)/%.o:	$(ROOTRIPDIR)/%.cc 
			rm -f $@ 
			@echo "Postponing $@ for separate Makefile" 

$(ROOTUTILLIB):	$(ROOTRIPOBJS)
		@cd rootstuff; $(MAKE) 

include $(SRCDIR)/Makefile 

FitManager.o:		FitManager.cc FitManager.hh wrkdir/ThrustFitManagerCUDA.o Variable.o 
			$(CXX) $(DEFINEFLAGS) $(CXXFLAGS) $(INCLUDES)  -o $@ -c $<

wrkdir/GooPdfCUDA.o:	wrkdir/CUDAglob.cu PdfBase.cu 
			$(CXX) $(CXXFLAGS) $(INCLUDES) -I. $(DEFINEFLAGS)  -c $< -o $@ 
			@echo "$@ done"

clean:
		@rm -f *.o wrkdir/*
		cd rootstuff; $(MAKE) clean 
