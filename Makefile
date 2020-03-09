CXX=g++
CPPFLAGS=-I$(MATHGEOLIB_INCLUDE_DIR) -I$(SOIL_INCLUDE_DIR)
CXXFLAGS=-Wall -ggdb
LDFLAGS=-L$(MATHGEOLIB_LIBRARY) -L$(SOIL_LIBRARY)
LIBS=-lMathGeoLib -lsoil -lGL -lGLU -lglut
EXE=ss

SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $(EXE)

clean:
	rm -f $(OBJECTS) $(EXE)