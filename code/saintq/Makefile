BOOST_INCLUDE_DIR = boost/include
DLIB_INCLUDE_DIR = dlib

CXXFLAGS += -O2 -std=gnu++0x -w

CPPFLAGS += -D NDEBUG
CPPFLAGS += -I $(DLIB_INCLUDE_DIR) -I $(BOOST_INCLUDE_DIR)

all: bin/saintq
bin/saintq: src/main.o
	$(LINK.cc) $^ $(LDLIBS) -o $@


.PHONY : clean
clean:
	$(RM) bin/saintq src/*.o
