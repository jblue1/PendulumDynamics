CXXFLAGS = -Wall  -I $(BOOST_ROOT) -std=c++17 -g -O2
Systems : systems.cpp systems.hpp
	$(CXX) $(CXXFLAGS) -c $<

Streaming_Observers : streaming_observers.cpp streaming_observer.hpp
	$(CXX) $(CXXFLAGS) -c $<

p_section_data.o : p_section_data.cpp systems.hpp streaming_observers.hpp
	h5c++ $(CXXFLAGS) -c $<

p_section_data : systems.o streaming_observers.o p_section_data.o
	h5c++ $(CXXFLAGS) -o ../bin/p_section_data $^

solve_system.o : solve_system.cpp systems.hpp streaming_observers.hpp
	$(CXX) $(CXXFLAGS) -c $<

solve_system : systems.o streaming_observers.o solve_system.o
	$(CXX) $(CXXFLAGS) -o ../bin/solve_system $^

clean:
	rm *.o
