# PTrees 
# Makefile
# Edoardo Carlesi 2018

include Makefile.config

EXEC=PTrees

PTrees :
	cd src/; make ${EXEC}
	mv src/${EXEC} bin/

clean  :
	cd src; rm *.o
