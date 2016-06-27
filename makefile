
ASSET_DIRS = $(shell find ./program_files/ -type d)
ASSET_FILES = $(shell find ./program_files/ -type f -name '*')

all: cnstools

cnstools: ./program_files $(ASSET_DIRS) $(ASSET_FILES)
	cd program_files; zip -r ../cnstools.zip *; cd ..
	echo '#!/usr/bin/env python' | cat - cnstools.zip > cnstools
	chmod +x cnstools
	rm cnstools.zip