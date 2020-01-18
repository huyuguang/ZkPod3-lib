#sudo apt-get install libomp-dev

CXX11_HOME = /usr

HOME = ../linux
THIRD_PARTY_HOME = ../thirdparty
YAS_HOME = ../thirdparty/yas/include
CSV_HOME = ../thirdparty/csv
PUB_HOME = ../public
SNARK_HOME = ../depends/libsnark
FF_HOME = $(SNARK_HOME)/depends/libff
MCL_HOME = $(SNARK_HOME)/depends/mcl/include

CXXFLAGS := \
 -g3 \
 -fPIC \
 -std=c++17 \
 -Wall \
 -Wextra \
 -gdwarf-2 \
 -gstrict-dwarf \
 -Wno-parentheses \
 -Wdeprecated-declarations \
 -fmerge-all-constants  \
 -march=native \
 -mtune=native \
 -I $(CXX11_HOME)/include \
 -isystem $(SNARK_HOME) \
 -isystem $(MCL_HOME) \
 -isystem $(FF_HOME) \
 -isystem $(YAS_HOME) \
 -isystem $(CSV_HOME) \
 -isystem $(THIRD_PARTY_HOME) \
 -I $(PUB_HOME) \
 -I. \
 -I.. \

MCL_FLAGS := \
 -DMCL_DONT_USE_OPENSSL

 BOOST_FLAGS := \
 -DBOOST_SPIRIT_THREADSAFE 

SNARK_FLAGS = \
 -DBN_SUPPORT_SNARK=1 \
 -DNO_PT_COMPRESSION=1 \
 -DCURVE_MCL_BN128 \
 -DMULTICORE=1 \
 -DUSE_ASM \
 -DNO_PROCPS \

RELEASE_FLAGS := \
 -O3 \
 -DNDEBUG \
 -DLINUX \

DEBUG_FLAGS := \
 -O0 \
 -D_DEBUG \
 -DLINUX \
 
LDFLAGS += \
 -fmerge-all-constants

LIBS := \
 -lpthread \
 -lrt \
 -lgmpxx \
 -lgmp \
 -lcryptopp \
 -lboost_filesystem \
 -lboost_system \
 -lboost_iostreams \
 -lboost_program_options \
 -ltbbmalloc_proxy \
 -ltbbmalloc \

ifdef USE_TBB
LIBS += -ltbb
CXXFLAGS += -DUSE_TBB
endif

BUILD := ./build/self
BUILDPUB := ./build/public

$(shell mkdir -p $(HOME)/bin)
$(shell mkdir -p $(BUILD))
$(shell mkdir -p $(BUILDPUB)/ecc)

DIR := . 

SRC := $(foreach d, $(DIR), $(wildcard $(d)/*.cc))

SRC += ../public/ecc/ecc.cc

RELEASE_OBJ := $(patsubst %.cc, %.o, $(SRC))

DEBUG_OBJ := $(patsubst %.cc, %.d.o, $(SRC))

RELEASE_BIN := pod_dummy
DEBUG_BIN := pod_dummyd

CXX := export LD_LIBRARY_PATH=${CXX11_HOME}/lib; ${CXX11_HOME}/bin/g++

#CC := export LD_LIBRARY_PATH=${CXX11_HOME}/lib; ${CXX11_HOME}/bin/gcc

all: $(RELEASE_BIN) $(DEBUG_BIN)

%.o : %.cc
	$(CXX) -c $^ $(CXXFLAGS) $(RELEASE_FLAGS) $(MCL_FLAGS) $(SNARK_FLAGS) $(BOOST_FLAGS) -o $(BUILD)/$@


%.d.o : %.cc
	$(CXX) -c $^ $(CXXFLAGS) $(DEBUG_FLAGS) $(MCL_FLAGS) $(SNARK_FLAGS) $(BOOST_FLAGS) -o $(BUILD)/$@



$(DEBUG_BIN) : $(DEBUG_OBJ)
	$(CXX) $(addprefix $(BUILD)/, $^) -o $(HOME)/bin/$@ $(LDFLAGS) $(LIBS)

$(RELEASE_BIN) : $(RELEASE_OBJ)
	$(CXX) $(addprefix $(BUILD)/, $^) -o $(HOME)/bin/$@ $(LDFLAGS) $(LIBS)

clean:
	find . -name "*.o" -delete
	rm -rf $(HOME)/bin/$(RELEASE_BIN) $(HOME)/bin/$(DEBUG_BIN)
	rm -rf build

.PHONY :
