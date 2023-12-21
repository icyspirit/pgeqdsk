#!/usr/bin/make

.PHONY: all clean

MODULENAME=fgeqdsk
SOURCE=fgeqdsk.f90

all: $(SOURCE)
	f2py -m $(MODULENAME) -c $(SOURCE)

clean:
	rm -f $(MODULENAME)$(shell python3-config --extension-suffix)
