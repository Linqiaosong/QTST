CXXFLAGS = -O3
CXXFLAGS += -std=c++11
RM = rm -f
CXX = g++
LD = g++
OBJS = main.o QDApplication.o QDynamic.o
.PHONY: clean
QTST : $(OBJS)
	$(LD) -static $(OBJS) -o QTST
main.o : main.cpp QDApplication.h QDynamic.h
	$(CXX) -c $(CXXFLAGS) main.cpp
QDApplication.o : QDApplication.cpp QDApplication.h QDynamic.h
	$(CXX) -c $(CXXFLAGS) QDApplication.cpp
QDynamic.o : QDynamic.cpp QDynamic.h
	$(CXX) -c $(CXXFLAGS) QDynamic.cpp
clean:
	$(RM) QTST $(OBJS)