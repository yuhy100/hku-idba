all:
	cd src/libs && make all
	cd src/tools && make all

clean:
	cd src/libs && make clean
	cd src/tools && make clean
	rm -f script/*~
	rm -rf idba-0.9

release:
	mkdir -p idba-0.9/src/tools
	mkdir -p idba-0.9/src/libs
	mkdir -p idba-0.9/ext-tool
	mkdir -p idba-0.9/bin
	mkdir -p idba-0.9/script
	cp Makefile idba-0.9
	cp src/tools/*.cpp idba-0.9/src/tools
	cp src/tools/Makefile idba-0.9/src/tools
	rm -f idba-0.9/umea.cpp
	cp src/libs/*.cpp idba-0.9/src/libs
	cp src/libs/*.h idba-0.9/src/libs
	cp src/libs/Makefile idba-0.9/src/libs
	cp ext-tool/* idba-0.9/ext-tool
	cp script/runIdba script/validate idba-0.9/script
	cp ReadMe.txt idba-0.9
	tar -czf idba-0.9.tar.gz idba-0.9