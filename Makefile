# 	METRO C++:	 MErgerTRees On C++ 
# A scalable C++/MPI merger tree code for cosmological simulations 
# Makefile
# Edoardo Carlesi 2018, 2019

include Makefile.config

EXEC_T=test
EXEC=MetroCPP

MetroCPP :
	cd src/; make ${EXEC}

# Rename the executable in zoom mode
ifeq ($(ZOOM_MODE), "true")
	mv bin/$(EXEC) bin/MetroCPP-Zoom
endif

test :
	cd src/; make ${EXEC_T}

.PHONY: all

all: ${EXEC} test


.PHONY: clean

clean  :
	cd src/; make clean;
