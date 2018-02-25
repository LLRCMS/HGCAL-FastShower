#makefile 


CC   =   g++
CCVERSIONGTEQ48 := $(shell expr `g++ -dumpversion | cut -f1,2 -d.` \>= 4.8)

#UCFLAGS = -O0 -g3 -Wall -gstabs+ -DSTANDALONE
UCFLAGS = -O3 -Wall -DSTANDALONE

# RUCFLAGS := $(shell root-config --cflags) -I./interface/ -I/usr/include/python2.6/ -std=c++14
# LIBS :=  $(shell root-config --libs) -lpython2.6 -lboost_python 


PYTHONPATH :=$(shell dirname `which python`)/../

RUCFLAGS := $(shell root-config --cflags) -I./interface/  -I/include/python2.7 -std=c++14
LIBS := $(shell root-config --libs) -L/opt/rh/python27/root/usr/lib64 -lpython2.7 -lboost_python 

vpath %.cpp ./src/json
vpath %.cpp ./src
vpath %.cpp ./bin

SRCPP = main.cpp\
	Cell.cpp\
	Point.cpp\
	Rectangle.cpp\
	Tree.cpp\
	Geometry.cpp\
	Generator.cpp\
	Parameters.cpp\
	ShowerParametrization.cpp\
	ShowerShape.cpp\
	ShowerShapeHexagon.cpp\
	ShowerShapeTriangle.cpp\
	Event.cpp\
	OutputService.cpp\
	jsoncpp.cpp\

#OBJCPP = $(SRCPP:.cpp=.o)
OBJCPP = $(patsubst %.cpp,lib/%.o,$(SRCPP))


ifeq "$(CCVERSIONGTEQ48)" "0"
  $(error Requires g++ version >= 4.8)
endif

all : bin/shower_simulation.exe

lib/%.o : %.cpp
	@echo "> compiling $*"
	@mkdir -p lib/
	@$(CC) -c $< $(UCFLAGS) $(RUCFLAGS) -o $@

bin/shower_simulation.exe : $(OBJCPP)
	@echo "> linking"
	$(CC) $^ $(LIBS) -o $@


clean:
	@echo "> Cleaning object files"
	@rm  -f lib/*.o
        
cleanall: clean
	@echo "> Cleaning executable"
	@rm -f bin/shower_simulation.exe
