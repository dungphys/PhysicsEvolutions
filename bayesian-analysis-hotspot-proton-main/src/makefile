CXX = g++
CXXFLAGS = -Wall -g
LDLIBS = -L/home/user/myDir -lgsl -lgslcblas -lm

RM=rm -f

BIN = /home/anh-dung/projects/current_projects/hot_spot_NR/bin
SRC = /home/anh-dung/projects/current_projects/hot_spot_NR/src

OBJS = root.o nintlib.o aux.o cross_sections.o
HPPS = root.hpp nintlib.hpp aux.hpp	cross_sections.hpp
SRCS = root.cpp nintlib.cpp aux.cpp	cross_sections.cpp	##root.cpp aux.cpp cross_sections.cpp

all: dsigma dsigma_pos test_rcut dsigma_test

test_rcut: test_rcut.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(BIN)/run_test test_rcut.o $(OBJS) $(LDLIBS)
test_rcut.o: test_rcut.cpp $(HPPS)
	$(CXX) $(CXXFLAGS) -c test_rcut.cpp

dsigma: dsigma.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(BIN)/run_prior dsigma.o $(OBJS) $(LDLIBS)
dsigma.o: dsigma.cpp $(HPPS)
	$(CXX) $(CXXFLAGS) -c dsigma.cpp

dsigma_pos: dsigma_pos.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(BIN)/run_pos dsigma_pos.o $(OBJS) $(LDLIBS)
dsigma_pos.o: dsigma_pos.cpp $(HPPS)
	$(CXX) $(CXXFLAGS) -c dsigma_pos.cpp
	
dsigma_test: dsigma_test.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(BIN)/run_test_dsigma dsigma_test.o $(OBJS) $(LDLIBS)
dsigma_test.o: dsigma_test.cpp $(HPPS)
	$(CXX) $(CXXFLAGS) -c dsigma_test.cpp

nintlib.o: nintlib.cpp nintlib.hpp
	$(CXX) $(CXXFLAGS) -c nintlib.cpp	
root.o: root.cpp root.hpp
	$(CXX) $(CXXFLAGS) -c root.cpp
aux.o: aux.cpp aux.hpp root.hpp nintlib.hpp
	$(CXX) $(CXXFLAGS) -c aux.cpp
cross_sections.o: cross_sections.cpp cross_sections.hpp aux.hpp root.hpp nintlib.hpp
	$(CXX) $(CXXFLAGS) -c cross_sections.cpp

	  
clean:
	$(RM) $(OBJS)
