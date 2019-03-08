
CXXFLAGS = -O3 -Wall -std=c++0x -I../nlopt/include -I../boost/include
all: SAINTexpress-int
SAINTexpress-int: main.o Stats.o BaitClass.o PreyClass.o InterClass.o UIClass.o Mapping.o
	$(CXX) $(CXXFLAGS) main.o Stats.o BaitClass.o PreyClass.o InterClass.o UIClass.o Mapping.o ../boost/lib/libboost_program_options.a -L../nlopt/lib -lnlopt -o ../bin/SAINTexpress-int

clean:
	$(RM) -rf *.o *~

