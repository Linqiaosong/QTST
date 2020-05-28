CXXFLAGS = -O3
CXXFLAGS += -std=c++11
CXXFLAGS += -fexec-charset=gbk
RM = del
CXX = g++
LD = g++
OBJS = main.o QDApplication.o QDynamic.o
.PHONY: clean
QTST.exe : $(OBJS)
	$(LD) -static $(OBJS) -o QTST.exe
main.o : main.cpp QDApplication.h QDynamic.h
	$(CXX) -c $(CXXFLAGS) main.cpp
QDApplication.o : QDApplication.cpp QDApplication.h QDynamic.h
	$(CXX) -c $(CXXFLAGS) QDApplication.cpp
QDynamic.o : QDynamic.cpp QDynamic.h
	$(CXX) -c $(CXXFLAGS) QDynamic.cpp
clean:
	$(RM) QTST.exe $(OBJS)