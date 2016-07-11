
PROGRAM_DIRECTORIES = $(shell find ./program_files/ -type d)
PROGRAM_FILES = $(shell find ./program_files/ -type f -name '*')

all: cnstools

cnstools: ./program_files $(PROGRAM_DIRECTORIES) $(PROGRAM_FILES)
	cd program_files; zip -r ../cnstools.zip *; cd ..
	echo '#!/usr/bin/env python' | cat - cnstools.zip > cnstools
	chmod +x cnstools
	rm cnstools.zip