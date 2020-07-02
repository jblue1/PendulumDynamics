Systems : systems.cpp systems.hpp
	g++ -c -O systems.cpp

Streaming_Observers : streaming_observers.cpp streaming_observer.cpp
	g++ -c -O streaming_observer.cpp

p_section_data.o : p_section_data.cpp systems.hpp streaming_observers.hpp
	h5c++ -std=c++17 -c -O p_section_data.cpp

p_section_data : systems.o streaming_observers.o p_section_data.o
	h5c++ -std=c++17 -O3 -o ../bin/p_section_data systems.o streaming_observers.o p_section_data.o

clean:
	rm *.o
