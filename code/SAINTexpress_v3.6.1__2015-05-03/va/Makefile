# test if std::begin and std::end exist, if they exist, va.hpp is an empty file.
CXXFLAGS = -std=c++0x
all: va.hpp
va.hpp:
	cat va_.hpp > va.hpp
	$(CXX) $(CXXFLAGS) v.cpp -o /dev/null && echo > va.hpp
