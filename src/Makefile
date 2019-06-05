# Makefile for METRO C++ 
# MErger TRees On C++
# Edoardo Carlesi 2018, 2019

include ../Makefile.config

# compiler optimization:
#CXXFLAGS += -g -O0           # debug mode
CXXFLAGS += -O3 -std=c++11 -fpermissive               # normal mode
#CXXFLAGS += -Wunused-result -Wsign-compare -Wunused-but-set-variable

EXEC=MetroCPP
EXEC_T=test

#=============================================================================#

SOURCE =\
	IOSettings.cpp Communication.cpp global_vars.cpp \
	Halo.cpp Grid.cpp MergerTree.cpp \
	Cosmology.cpp utils.cpp	spline.cpp

OBJS  =  $(SOURCE:.cpp=.o)

# clear out all suffixes
.SUFFIXES:
# list only those we use
.SUFFIXES: .o .cpp

#=============================================================================#

.cpp.o :
	$(CXX) $(CXXFLAGS) $(OPT) -c $<

.PHONY: all clean

all : 	$(EXEC)

MetroCPP : $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPT) -o $(EXEC) main.cpp $(OBJS) $(LDFLAGS) $(LDLIBS)
	mv $(EXEC) ../bin/

# Check if file format *.o exists before cleaning
ifneq ("$(wildcard *.o)","")
FILE_EXISTS = 1
else
FILE_EXISTS = 0
endif

clean :
	@if [ -f $(EXEC) ]; then rm $(EXEC); fi;
	@if [ -f $(EXEC_T) ]; then rm $(EXEC_T); fi;
	@if [ -f ../bin/$(EXEC) ]; then rm ../bin/$(EXEC); fi;
	@if [ -f ../bin/$(EXEC_T) ]; then rm ../bin/$(EXEC_T); fi;
	@if [ "$(FILE_EXISTS)" -eq "1" ]; then rm *.o; fi;
