# "make print-X print-XS print-Y" prints values of macros X, XS and Y
# see http://www.drdobbs.com/tools/debugging-makefiles/197003338
print-%: ; @echo $* is $($*)

HOST_SYSTEM = $(shell uname | cut -f 1 -d_)
SYSTEM ?= $(HOST_SYSTEM)

NAME = sensitivity

# Use PSUFFIX for executable, for debugging under XCode, cannot have suffix,
# so leave PSUFFIX blank for debugging

ifeq ($(SYSTEM),Darwin)
SUFFIX  = .osx
PSUFFIX =
endif

ifeq ($(SYSTEM),Linux)
SUFFIX  = .il
PSUFFIX = .il
endif

PROG    = $(NAME)$(PSUFFIX)
DEPEND  = src/dependencies$(SUFFIX)

CXXFILES   =  $(NAME).cc Individual.cc Population.cc SumStat.cc Performance.cc
OBJFILES   = $(CXXFILES:.cc=.o)

# defs for linking to sim_client.cc instead of main-alone.cc
# src/sim_client.cc is symbolic link to original in grpc directory
# proto directory
OCLIENT    = sim_client.o
CXXCLIENT  = src/sim_client.cc
CLIENT     = $(NAME)_client$(SUFFIX)
DEFS       = -DCLIENT_LINK
PROTO_PATH = proto
POBJFILES = $(PROTO_PATH)/simcontrol.pb.o $(PROTO_PATH)/simcontrol.grpc.pb.o
PLDFLAGS = `pkg-config --libs protobuf grpc++ grpc` -lgrpc++_reflection

APPHEAD   = $(NAME).h

# Google performance tools, see https://github.com/gperftools/gperftools
# add -lprofiler -ltcmalloc\ to LDFLAGS
# mkdir profile
# turn on with environment variables
# setenv HEAPCHECK normal; setenv HEAPPROFILE profile/heap.prof; setenv CPUPROFILE profile/cpu.prof
# get output with pprof, e.g., pprof --text --cum sensitivity.osx profile/*
# or pprof, e.g., pprof --pdf --cum sensitivity.osx profile/* > prof.pdf

# find misaligned bugs: -fsanitize=undefined -fno-omit-frame-pointer in compile

CXX = g++
CXXFLAGS += -std=c++17      # c++1z also, but can now use c++17 explicitly
CXXFLAGS += -pedantic -Wall -Wextra -W \
  -Wconversion -Wshadow -Wundef \
  -Wpointer-arith -Wcast-align \
  -Wwrite-strings -Wstrict-prototypes \
  -Wcast-qual -Wconversion \
  -g $(INCFLAGS) $(DEFS) -O3 #-pg #-DDEBUG
LDFLAGS += -L$(HOME)/sim/simlib/lib_osx -L/opt/local/lib -lfmt\
             -lutilSAF -lboost_system-mt -lboost_filesystem-mt -lgsl -lgslcblas\

INCFLAGS = -I$(HOME)/sim/simlib/include -Isrc -I$(PROTO_PATH) -isystem /opt/local/include
VPATH 	= src:$(PROTO_PATH)
GARBAGE = core error.* 
PWD := $(shell pwd)

.PHONY: $(NAME).first
$(NAME).first:
	$(MAKE) $(PROG) $(MFLAGS) "CXXFLAGS = $(CXXFLAGS) -DAPPL_H=\\\"$(APPHEAD)\\\""

$(PROG): $(OBJFILES) main-alone.o
	$(CXX) $(CXXFLAGS) -o $@ $(OBJFILES) main-alone.o $(LDFLAGS)

debug: $(PROG)
	dsymutil $(PROG)

client: proto $(CLIENT)

$(CLIENT): $(OBJFILES) $(OCLIENT)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJFILES) $(OCLIENT) $(POBJFILES) $(LDFLAGS) $(PLDFLAGS)

.PHONY: proto
proto:
	cd $(HOME)/sim/grpcControl; $(MAKE) proto

depend:
	gcc -MM $(CXXFLAGS)  -DAPPL_H=\"$(APPHEAD)\" $(addprefix src/, $(CXXFILES)) src/main-alone.cc $(CXXCLIENT)> $(DEPEND)
    #perl -p -i -e 's/^(\S)/src\/\1/' src/dependencies.osx   # prepend 'src/' for targets

# must update dependency file by typing "make depend"

include $(DEPEND)

clean:
	-rm -f  *.o $(PROG) $(CLIENT) $(GARBAGE)
	-rm -rf *.dSYM

cleanproto: 
	cd $(HOME)/sim/grpcControl; $(MAKE) cleanobj
	
cleanall: clean cleanproto

help:
	@echo '  make $(NAME) -  to make the application for running'
	@echo '  make clean -    to remove all files but the source'
	@echo '  reset flags for debugging or profiling'
