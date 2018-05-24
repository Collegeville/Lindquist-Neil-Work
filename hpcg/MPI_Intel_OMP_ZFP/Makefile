# -*- Makefile -*-

arch = MPI_Intel_OMP_ZFP
setup_file = setup/Make.$(arch)

include $(setup_file)


HPCG_DEPS = src/GenerateGeometry.o \
	    src/MixedBaseCounter.o \
	    src/ReadHpcgDat.o \
	    src/YAML_Doc.o \
	    src/YAML_Element.o \
	    src/mytimer.o \
	    src/ComputeOptimalShapeXYZ.o \
	    src/CheckAspectRatio.o \
	    src/OutputFile.o \
	    src/init.o \
	    src/finalize.o

# These header files are included in many source files, so we recompile every file if one or more of these header is modified.
PRIMARY_HEADERS = $(INCdir)/Geometry.hpp $(INCdir)/SparseMatrix.hpp $(INCdir)/Vector.hpp \
                  $(INCdir)/CGData.hpp $(INCdir)/MGData.hpp $(INCdir)/hpcg.hpp

all: bin/xhpcg

bin/xhpcg: src/main.o $(HPCG_DEPS)
	$(LINKER) $(LINKFLAGS) src/main.o $(HPCG_DEPS) $(HPCG_LIBS) -o bin/xhpcg

clean:
	rm -f src/*.o bin/xhpcg

.PHONY: all clean

src/main.o: $(SRCdir)/main.cpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/GenerateGeometry.o: $(SRCdir)/GenerateGeometry.cpp $(SRCdir)/GenerateGeometry.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/MixedBaseCounter.o: $(SRCdir)/MixedBaseCounter.cpp $(SRCdir)/MixedBaseCounter.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ReadHpcgDat.o: $(SRCdir)/ReadHpcgDat.cpp $(SRCdir)/ReadHpcgDat.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/YAML_Doc.o: $(SRCdir)/YAML_Doc.cpp $(SRCdir)/YAML_Doc.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/YAML_Element.o: $(SRCdir)/YAML_Element.cpp $(SRCdir)/YAML_Element.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/finalize.o: $(SRCdir)/finalize.cpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/init.o: $(SRCdir)/init.cpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/mytimer.o: $(SRCdir)/mytimer.cpp $(SRCdir)/mytimer.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeOptimalShapeXYZ.o: $(SRCdir)/ComputeOptimalShapeXYZ.cpp $(SRCdir)/ComputeOptimalShapeXYZ.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/CheckAspectRatio.o: $(SRCdir)/CheckAspectRatio.cpp $(SRCdir)/CheckAspectRatio.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/OutputFile.o: $(SRCdir)/OutputFile.cpp $(SRCdir)/OutputFile.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

