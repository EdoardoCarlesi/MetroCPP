# 	METRO C++:	 MErgerTRees On C++ 
# A scalable C++/MPI merger tree code for cosmological simulations 
# Makefile
# Edoardo Carlesi 2018

include Makefile.config

EXEC=MetroCPP
EXEC_T=test

MetroCPP :
	cd src/; make ${EXEC}

test :
	cd src/; make ${EXEC_T}

clean  :
	rm bin/*
	cd src; rm *.o; 
