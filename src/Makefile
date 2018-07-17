# Makefile for METRO C++ 
# MErger TRees On C++
# Edoardo Carlesi 2018

include ../Makefile.config

# compiler optimization:
#CXXFLAGS += -g -O0           # debug mode
CXXFLAGS += -O3 -std=c++11               # normal mode
CXXFLAGS += -Wunused-result -fpermissive
EXEC = MetroCPP
EXEC_T = test

#=============================================================================#

SOURCE =\
	IOSettings.cpp Communication.cpp global_vars.cpp \
	Halo.cpp Grid.cpp MergerTree.cpp \
	Gravity.cpp Cosmology.cpp utils.cpp

OBJS  =  $(SOURCE:.cpp=.o)

# clear out all suffixes
.SUFFIXES:
# list only those we use
.SUFFIXES: .o .cpp

#=============================================================================#

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: all clean

all : 	MetroCPP test

MetroCPP : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) main.cpp $(OBJS) $(LDFLAGS) $(LDLIBS)
	mv $(EXEC) ../bin/

test : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC_T) test.cpp $(OBJS) $(LDFLAGS) $(LDLIBS)
	mv $(EXEC_T) ../bin/

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
