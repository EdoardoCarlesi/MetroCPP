# Makefile
# Edoardo Carlesi 2018

include ../Makefile.config

# compiler optimization:
#CXXFLAGS += -g -O0           # debug mode
CXXFLAGS += -O3 -std=c++11               # normal mode
CXXFLAGS += -Wunused-result -fpermissive
EXEC = PTrees
EXEC_T = test

#=============================================================================#

SOURCE =\
	IOSettings.cpp Comm.cpp general.cpp \
	Particle.cpp Halo.cpp Grid.cpp MergerTree.cpp

OBJS  =  $(SOURCE:.cpp=.o)

# clear out all suffixes
.SUFFIXES:
# list only those we use
.SUFFIXES: .o .cpp

#=============================================================================#

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

#all :
#	make PTrees; make test

PTrees : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) main.cpp $(OBJS) $(LDFLAGS) $(LDLIBS)
	mv PTrees ../bin/

#test : $(OBJS)
#	$(CXX) $(CXXFLAGS) -o $(EXEC_T) test.cpp $(OBJS) $(LDFLAGS) $(LDLIBS)


clean :
	rm *.o;	rm $(EXEC); rm $(EXEC_T) 
	rm ../bin/$(EXEC); rm ../bin/$(EXEC_T) 
