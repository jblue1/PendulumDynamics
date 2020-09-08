CXXFLAGS = -Wall  -I /Users/johnblue/Programs/Libraries/boost_1_73_0 -std=c++17 -g
LDFLAGS  = -L /Users/johnblue/Programs/Libraries/boost_1_73_0/stage/lib/ -rpath /Users/johnblue/Programs/Libraries/boost_1_73_0/stage/lib/
LDLIBS   = -l boost_unit_test_framework

systems : ../src/systems.cpp ../src/systems.hpp
	$(CXX) $(CXXFLAGS) -c -O ../src/systems.cpp

systems-test.o: systems-test.cpp ../src/systems.hpp
	$(CXX) $(CXXFLAGS) -c $<

systems-test: ../src/systems.o systems-test.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -o $@ $^
