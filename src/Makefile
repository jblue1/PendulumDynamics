CXXFLAGS = -Wall  -I /Users/johnblue/Programs/Libraries/boost_1_73_0 -std=c++17 -g -O2
Systems : systems.cpp systems.hpp
	$(CXX) $(CXXFLAGS) -c $<

Streaming_Observers : streaming_observers.cpp streaming_observer.hpp
	$(CXX) $(CXXFLAGS) -c $<

p_section_data.o : p_section_data.cpp systems.hpp streaming_observers.hpp
	h5c++ $(CXXFLAGS) -c $<

p_section_data : systems.o streaming_observers.o p_section_data.o
	h5c++ $(CXXFLAGS) -o ../bin/p_section_data $^

solve_system.o : solve_system.cpp systems.hpp streaming_observers.hpp
	g++ -std=c++17 -c -O solve_system.cpp

solve_system : systems.o streaming_observers.o solve_system.o
	g++ -std=c++17 -O3 -o ../bin/solve_system systems.o streaming_observers.o solve_system.o

clean:
	rm *.o
