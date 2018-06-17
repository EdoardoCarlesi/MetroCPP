# Makefile
# Edoardo Carlesi 2018

include ../Makefile.config

# compiler optimization:
#CXXFLAGS += -g -O0           # debug mode
CXXFLAGS += -O3               # normal mode
EXEC = PTrees

#=============================================================================#

SOURCE = Particle.cpp Halo.cpp Grid.cpp \
	IO.cpp general.cpp

OBJS  =  $(SOURCE:.cpp=.o)

# clear out all suffixes
.SUFFIXES:
# list only those we use
.SUFFIXES: .o .cpp

#=============================================================================#

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

PTrees : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) main.cpp $(OBJS) $(LDFLAGS) $(LDLIBS)

clean :
	rm *.o
