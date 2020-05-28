CXXFLAGS = -O3
CXXFLAGS += -std=c++11
CXXFLAGS += -static
RM = rm -f
CXX = clang++
LD = clang++
OBJS = main.o QDApplication.o QDynamic.o
.PHONY: clean
QTST : $(OBJS)
	$(LD) $(OBJS) -o QTST
main.o : main.cpp QDApplication.h QDynamic.h
	$(CXX) -c $(CXXFLAGS) main.cpp
QDApplication.o : QDApplication.cpp QDApplication.h QDynamic.h
	$(CXX) -c $(CXXFLAGS) QDApplication.cpp
QDynamic.o : QDynamic.cpp QDynamic.h
	$(CXX) -c $(CXXFLAGS) QDynamic.cpp
clean:
	$(RM) QTST $(OBJS)