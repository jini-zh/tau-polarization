CXXFLAGS ?= -O2 -pipe

rootcmp: rootcmp.cpp
	g++ $(CXXFLAGS) `root-config --cflags --ldflags --libs` $< -o $@
