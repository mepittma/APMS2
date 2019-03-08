all: bin/SAINTexpress-spc bin/SAINTexpress-int
bin/SAINTexpress-int: boost nlopt va/va.hpp
	cd SAINT-MRF-int/; $(MAKE)
bin/SAINTexpress-spc: boost nlopt va/va.hpp
	cd SAINT-MRF-spc/; $(MAKE)
va/va.hpp:
	cd va/; $(MAKE) -i
nlopt:
	cd nloptsrc/; ./configure --prefix=$(CURDIR)/nlopt; $(MAKE) install
boost: 
	mkdir -p boost/include; ln -s ../../boostsrc/boost boost/include/boost; cd boostsrc/; ./bootstrap.sh; ./b2 -j 4 --with-program_options link=static stage --stagedir=../boost/
