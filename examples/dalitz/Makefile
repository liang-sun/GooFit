#------------------------------------------------------------------------------
CXX=nvcc
LD=g++  
OutPutOpt = -o

#CXXFLAGS     = -w -std=c++11 -O3 #-arch=sm_20 -g -G 
CXXFLAGS     = -w -std=c++11 -D_MWAITXINTRIN_H_INCLUDED -D_FORCE_INLINES -O3 #-g -G 
ifeq ($(CUDALOCATION), )
CUDALOCATION = /usr/local/cuda/
endif
CUDAHEADERS = $(CUDALOCATION)/include/

GOOFITDIR = ${PWD}/../../
#EIGENINCDIR=$(HOME)/MinuitThrust/eigen3
SRCDIR = $(GOOFITDIR)/PDFs
WRKDIR = $(GOOFITDIR)/wrkdir/

ifneq ($(TARGET_OMP),)
CXXFLAGS    += -Xcompiler -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP -lgomp
LDFLAGS     += -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP -lgomp
endif 

INCLUDES += -I$(CUDAHEADERS) -I$(SRCDIR) -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff 
#INCLUDES += -I$(EIGENINCDIR) -I$(GOOFITDIR)/MCBooster
INCLUDES += -I$(GOOFITDIR)/MCBooster
LIBS += -L$(CUDALOCATION)/lib64 -lcudart -L$(GOOFITDIR)/rootstuff -lRootUtils 
LIBS += -L$(GOOFITDIR) -lGooFit 

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ROOT_INCLUDES = -I/usr/include/root
ROOT_LIBS     =  -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lMatrix -lPhysics -lMathCore -pthread -lThread -lMinuit -rdynamic -lFoam 

#THRUSTO		= $(WRKDIR)/Variable.o $(WRKDIR)/FitManager.o $(WRKDIR)/GooPdfCUDA.o $(WRKDIR)/Faddeeva.o $(WRKDIR)/FitControl.o $(WRKDIR)/PdfBase.o $(WRKDIR)/DataSet.o $(WRKDIR)/BinnedDataSet.o $(WRKDIR)/UnbinnedDataSet.o $(WRKDIR)/FunctorWriter.o 
GOOFITLIB	= $(GOOFITDIR)/libGooFit.a
ROOTRIPDIR	= $(GOOFITDIR)/rootstuff
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 
LINKOBJ         = link.o

.SUFFIXES: 

all:	dalitz

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -dc -o $@ $<

dalitz:		DalitzFit.o  $(GOOFITLIB) $(ROOTUTILLIB) 
		$(CXX) -dlink $^ -o $(LINKOBJ)
		$(LD) $(LDFLAGS) $^ $(LINKOBJ) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz_plot:		DalitzFit_plot.o  $(GOOFITLIB) $(ROOTUTILLIB) 
		$(CXX) -dlink $^ -o $(LINKOBJ)
		$(LD) $(LDFLAGS) $^ $(LINKOBJ) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz2:		DalitzFit2.o $(GOOFITLIB) $(ROOTUTILLIB) 
		$(CXX) -dlink $^ -o $(LINKOBJ)
		$(LD) $(LDFLAGS) $^ $(LINKOBJ) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz2_gen:		DalitzFit2_gen.o $(GOOFITLIB) $(ROOTUTILLIB) 
		$(CXX) -dlink $^ -o $(LINKOBJ)
		$(LD) $(LDFLAGS) $^ $(LINKOBJ) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz3:		DalitzFit3.o $(GOOFITLIB) $(ROOTUTILLIB) 
		$(CXX) -dlink $^ -o $(LINKOBJ)
		$(LD) $(LDFLAGS) $^ $(LINKOBJ) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz_symm:		DalitzFit_symm.o $(GOOFITLIB) $(ROOTUTILLIB) 
		$(CXX) -dlink $^ -o $(LINKOBJ)
		$(LD) $(LDFLAGS) $^ $(LINKOBJ) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz_D3pi_sig:		DalitzFit_D3pi_sig.o $(GOOFITLIB) $(ROOTUTILLIB) 
		$(CXX) -dlink $^ -o $(LINKOBJ)
		$(LD) $(LDFLAGS) $^ $(LINKOBJ) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz_Ds3pi_sig:		DalitzFit_Ds3pi_sig.o $(THRUSTO) $(ROOTUTILLIB) 
		$(CXX) -dlink $^ -o $(LINKOBJ)
		$(LD) $(LDFLAGS) $^ $(LINKOBJ) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

dalitz_D3K_sig:		DalitzFit_D3K_sig.o $(GOOFITLIB) $(ROOTUTILLIB) 
		$(CXX) -dlink $^ -o $(LINKOBJ)
		$(LD) $(LDFLAGS) $^ $(LINKOBJ) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
		@echo "$@ done"

clean:
		@rm -f *.o core dalitz dalitz.o 
