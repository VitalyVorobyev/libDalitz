Sources=EvtComplex.cpp EvtConst.cpp EvtKine.cpp EvtResonance2.cpp EvtVector3C.cpp EvtVector3R.cpp EvtVector4C.cpp EvtVector4R.cpp EvtTensor4C.cpp symdalitzmodel.cpp dalitzmodel.cpp kspipimodel.cpp dalitzphasespace.cpp b0tod0pipimodel.cpp modelintegral.cpp libdalitz.cpp randomdalitzpoint.cpp dalitzmcintegral.cpp dalitzgenerator.cpp
Executable=libdalitz.so
CFlags=-c -Wall -fPIC -g -Iinc -std=c++11 -I. `root-config --cflags`
LDFlags= -shared -fPIC -std=c++11 -I. -Wl,--no-as-needed `root-config --glibs` -lm -lstdc++ -lRooFit -lRooFitCore -ldalitz
ObjectDir=obj/
SourceDir=src/
BinDir=bin/

CC=g++
RM=rm

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
Objects=$(Sources:.cpp=.o)
CSources=$(addprefix $(SourceDir),$(Sources))
CObjects=$(addprefix $(ObjectDir),$(Objects))
CExecutable=$(addprefix $(BinDir),$(Executable))

all: $(CSources) $(CExecutable)

$(CExecutable): $(CObjects)
	$(CC) $(LDFlags) $(CObjects) -o $@

$(ObjectDir)%.o: $(SourceDir)%.cpp
	$(CC) $(CFlags) $< -o $@

clean:
	$(RM) $(CObjects)
