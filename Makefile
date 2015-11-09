all: libdalitz.so

libdalitz.so: EvtComplex.o EvtConst.o EvtKine.o EvtResonance2.o EvtVector3C.o EvtVector3R.o EvtVector4C.o EvtVector4R.o EvtTensor4C.o symdalitzmodel.o dalitzmodel.o kspipimodel.o dalitzphasespace.o b0tod0pipimodel.o modelintegral.o libdalitz.o
	g++ -Wall -shared -fPIC -std=c++11 -o libdalitz.so EvtComplex.o EvtConst.o EvtKine.o EvtResonance2.o EvtVector3C.o EvtVector3R.o EvtVector4C.o EvtVector4R.o EvtTensor4C.o symdalitzmodel.o dalitzmodel.o kspipimodel.o dalitzphasespace.o b0tod0pipimodel.o modelintegral.o libdalitz.o -I. -lm -lstdc++ 

EvtComplex.o: EvtComplex.cpp
	g++ -fPIC -c EvtComplex.cpp

EvtConst.o: EvtConst.cpp
	g++ -fPIC -c EvtConst.cpp

EvtKine.o: EvtKine.cpp
	g++ -fPIC -c EvtKine.cpp

EvtResonance2.o: EvtResonance2.cpp
	g++ -fPIC -c EvtResonance2.cpp

EvtVector3C.o: EvtVector3C.cpp
	g++ -fPIC -c EvtVector3C.cpp

EvtVector3R.o: EvtVector3R.cpp
	g++ -fPIC -c EvtVector3R.cpp

EvtVector4C.o: EvtVector4C.cpp
	g++ -fPIC -c EvtVector4C.cpp

EvtVector4R.o: EvtVector4R.cpp
	g++ -fPIC -c EvtVector4R.cpp

EvtTensor4C.o: EvtTensor4C.cpp
	g++ -fPIC -c EvtTensor4C.cpp

dalitzmodel.o: dalitzmodel.cpp
	g++ -fPIC -c dalitzmodel.cpp

symdalitzmodel.o: symdalitzmodel.cpp
	g++ -fPIC -c symdalitzmodel.cpp

dalitzphasespace.o: dalitzphasespace.cpp
	g++ -fPIC -c dalitzphasespace.cpp

kspipimodel.o: kspipimodel.cpp
	g++ -fPIC -std=c++11 -c kspipimodel.cpp

b0tod0pipimodel.o: b0tod0pipimodel.cpp
	g++ -fPIC -std=c++11 -c b0tod0pipimodel.cpp

modelintegral.o: modelintegral.cpp
	g++ -fPIC -c modelintegral.cpp

libdalitz.o: libdalitz.cpp
	g++ -fPIC -c libdalitz.cpp

clean:
	rm -rf *.o libdalitz.so

