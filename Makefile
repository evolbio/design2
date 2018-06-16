# "make print-X print-XS print-Y" prints values of macros X, XS and Y
# see http://www.drdobbs.com/tools/debugging-makefiles/197003338
print-%: ; @echo $* is $($*)

HOST_SYSTEM = $(shell uname | cut -f 1 -d_)
SYSTEM ?= $(HOST_SYSTEM)

NAME = sensitivity

ifeq ($(SYSTEM),Darwin)
SUFFIX  = osx
endif

ifeq ($(SYSTEM),Linux)
SUFFIX  = il
endif

PROG    = $(NAME).$(SUFFIX)
DEPEND  = src/dependencies.$(SUFFIX)

CXXFILES   =  $(NAME).cc Individual.cc Population.cc SumStat.cc
OBJFILES   = $(CXXFILES:.cc=.o)

# defs for linking to sim_client.cc instead of main-alone.cc
# src/sim_client.cc is symbolic link to original in grpc directory
# proto directory
OCLIENT    = sim_client.o
CXXCLIENT  = src/sim_client.cc
CLIENT     = $(NAME)_client.$(SUFFIX)
DEFS       = -DCLIENT_LINK
PROTO_PATH = proto
POBJFILES = $(PROTO_PATH)/simcontrol.pb.o $(PROTO_PATH)/simcontrol.grpc.pb.o
PLDFLAGS = `pkg-config --libs protobuf grpc++ grpc` -lgrpc++_reflection

APPHEAD   = $(NAME).h

CXX = g++
CXXFLAGS += -std=c++17      # c++1z also, but can now use c++17 explicitly
CXXFLAGS += -pedantic -Wall -Wextra -W \
  -Wconversion -Wshadow -Wundef \
  -Wpointer-arith -Wcast-align \
  -Wwrite-strings -Wstrict-prototypes \
  -Wcast-qual -Wconversion \
  -g $(INCFLAGS) $(DEFS) -O3 #-pg #-DDEBUG
LDFLAGS += -L$(HOME)/sim/simlib/lib_osx -L/opt/local/lib\
             -lutilSAF -lboost_system-mt -lboost_filesystem-mt -lgsl -lgslcblas

INCFLAGS = -I$(HOME)/sim/simlib/include -Isrc -I$(PROTO_PATH) -isystem /opt/local/include
VPATH 	= src:$(PROTO_PATH)
GARBAGE = core error.* 
PWD := $(shell pwd)

.PHONY: $(NAME)
$(NAME):
	$(MAKE) $(PROG) $(MFLAGS) "CXXFLAGS = $(CXXFLAGS) -DAPPL_H=\\\"$(APPHEAD)\\\""

$(PROG): $(OBJFILES) main-alone.o
	$(CXX) $(CXXFLAGS) -o $@ $(OBJFILES) main-alone.o $(LDFLAGS)

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
	cd $(HOME)/sim/grpcControl; $(MAKE) cleanobj

help: 
	@echo '  make $(NAME) -  to make the application for running'
	@echo '  make clean -    to remove all files but the source'
	@echo '  reset flags for debugging or profiling'
