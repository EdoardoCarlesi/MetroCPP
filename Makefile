# PTrees 
# Makefile
# Edoardo Carlesi 2018

include Makefile.config

EXEC=PTrees
EXEC_T=test

PTrees :
	cd src/; make ${EXEC}

test :
	cd src/; make ${EXEC_T}

clean  :
	rm bin/*
	cd src; rm *.o; 
