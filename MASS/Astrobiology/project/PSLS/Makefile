
SOURCES=  pssl.py sls.py transit.py universal_pattern.py FortranIO.py 

all: build install

build: setup.py 
	python setup.py build	

install: setup.py 
	python setup.py install 
clean:
	python setup.py clean
	rm -fr build
