CXXFLAGS = -Wall  -I $(BOOST_ROOT) -std=c++17 -g -O2
Systems : systems.cpp systems.hpp
	$(CXX) $(CXXFLAGS) -c $<

p_section_data.o : p_section_data.cpp systems.hpp
	h5c++ $(CXXFLAGS) -c $<

p_section_data : systems.o p_section_data.o
	h5c++ $(CXXFLAGS) -o ../bin/p_section_data $^

solve_system.o : solve_system.cpp systems.hpp
	$(CXX) $(CXXFLAGS) -c $<

solve_system : systems.o solve_system.o
	$(CXX) $(CXXFLAGS) -o ../bin/solve_system $^

basins.o : basins.cpp systems.hpp
	h5c++ $(CXXFLAGS) -c $<

basins : systems.o basins.o
	h5c++ $(CXXFLAGS) -o ../bin/basins $^

clean:
	rm *.o
