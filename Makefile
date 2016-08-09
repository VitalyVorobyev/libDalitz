SourceDir=src/
ObjectDir=obj/
BinDir=bin/
Sources=$(notdir $(wildcard $(SourceDir)*.cpp))
Executable=libdalitz.so
CFlags=-c -Wall -fPIC -g -Iinc -std=c++14 
LDFlags= -shared -fPIC -std=c++14 -I. -lm -lstdc++
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
