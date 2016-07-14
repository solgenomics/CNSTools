
SOURCE_DIRECTORIES = $(shell find ./source/ -type d)
SOURCE_FILES = $(shell find ./source/ -type f -name '*')
DOC_DIRECTORIES = $(shell find ./docs/ -type d)
DOC_FILES = $(shell find ./docs/ -type f -name '*')
DOC_PATH_CHANGE_COMMAND = import os; import sys; sys.path.insert(0, os.path.abspath(\"../source\"));

all: cnstools docs

cnstools: macwash_source ./source $(SOURCE_DIRECTORIES) $(SOURCE_FILES)
	mkdir -p ./build/
	find ./build/ -name "cnstools" -delete
	cp -r ./source/cnstools ./build/cnstools
	find ./build/ -name "*.pyc" -delete
	cd ./build/cnstools/; zip -r ../cnstools.zip *; cd ../..
	rm -r ./build/cnstools
	echo '#!/usr/bin/env python' | cat - ./build/cnstools.zip > ./build/cnstools
	chmod +x ./build/cnstools
	rm ./build/cnstools.zip

macwash_source: ./source $(SOURCE_DIRECTORIES) $(SOURCE_FILES)
	find ./source/ -name "._*" -delete
	find ./source/ -name ".DS_Store" -delete

docs: macwash_source $(DOC_DIRECTORIES) $(DOC_FILES) $(SOURCE_DIRECTORIES) $(SOURCE_FILES)
	mkdir -p ./docs
	sphinx-apidoc -T -M -F -e -o ./docs ./source -H CNStools
	sed -i "/use\ os.path.abspath/c\$(DOC_PATH_CHANGE_COMMAND)" ./docs/conf.py
	cd docs; make html; cd ..