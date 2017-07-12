# -*- Makefile -*-

arch = MPI_gprof_OMP_Double
setup_file = setup/Make.$(arch)

include $(setup_file)


HPCG_DEPS = src/CG.o \
	    src/CG_ref.o \
	    src/TestCG.o \
	    src/ComputeResidual.o \
	    src/ExchangeHalo.o \
	    src/GenerateGeometry.o \
	    src/GenerateProblem.o \
	    src/GenerateProblem_ref.o \
	    src/CheckProblem.o \
	    src/MixedBaseCounter.o \
	    src/OptimizeProblem.o \
	    src/ReadHpcgDat.o \
	    src/ReportResults.o \
	    src/SetupHalo.o \
	    src/SetupHalo_ref.o \
	    src/TestSymmetry.o \
	    src/TestNorms.o \
	    src/WriteProblem.o \
	    src/YAML_Doc.o \
	    src/YAML_Element.o \
	    src/ComputeDotProduct.o \
	    src/ComputeDotProduct_ref.o \
	    src/mytimer.o \
	    src/ComputeOptimalShapeXYZ.o \
	    src/ComputeSPMV.o \
	    src/ComputeSPMV_ref.o \
	    src/ComputeSYMGS.o \
	    src/ComputeSYMGS_ref.o \
	    src/ComputeWAXPBY.o \
	    src/ComputeWAXPBY_ref.o \
	    src/ComputeMG_ref.o \
	    src/ComputeMG.o \
	    src/ComputeProlongation_ref.o \
	    src/ComputeRestriction_ref.o \
	    src/CheckAspectRatio.o \
	    src/OutputFile.o \
	    src/GenerateCoarseProblem.o \
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

src/CG.o: $(SRCdir)/CG.cpp $(SRCdir)/CG.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/CG_ref.o: $(SRCdir)/CG_ref.cpp $(SRCdir)/CG_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/TestCG.o: $(SRCdir)/TestCG.cpp $(SRCdir)/TestCG.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeResidual.o: $(SRCdir)/ComputeResidual.cpp $(SRCdir)/ComputeResidual.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ExchangeHalo.o: $(SRCdir)/ExchangeHalo.cpp $(SRCdir)/ExchangeHalo.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/GenerateGeometry.o: $(SRCdir)/GenerateGeometry.cpp $(SRCdir)/GenerateGeometry.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/GenerateProblem.o: $(SRCdir)/GenerateProblem.cpp $(SRCdir)/GenerateProblem.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/GenerateProblem_ref.o: $(SRCdir)/GenerateProblem_ref.cpp $(SRCdir)/GenerateProblem_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/CheckProblem.o: $(SRCdir)/CheckProblem.cpp $(SRCdir)/CheckProblem.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/MixedBaseCounter.o: $(SRCdir)/MixedBaseCounter.cpp $(SRCdir)/MixedBaseCounter.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/OptimizeProblem.o: $(SRCdir)/OptimizeProblem.cpp $(SRCdir)/OptimizeProblem.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ReadHpcgDat.o: $(SRCdir)/ReadHpcgDat.cpp $(SRCdir)/ReadHpcgDat.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ReportResults.o: $(SRCdir)/ReportResults.cpp $(SRCdir)/ReportResults.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/SetupHalo.o: $(SRCdir)/SetupHalo.cpp $(SRCdir)/SetupHalo.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/SetupHalo_ref.o: $(SRCdir)/SetupHalo_ref.cpp $(SRCdir)/SetupHalo_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/TestSymmetry.o: $(SRCdir)/TestSymmetry.cpp $(SRCdir)/TestSymmetry.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/TestNorms.o: $(SRCdir)/TestNorms.cpp $(SRCdir)/TestNorms.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/WriteProblem.o: $(SRCdir)/WriteProblem.cpp $(SRCdir)/WriteProblem.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/YAML_Doc.o: $(SRCdir)/YAML_Doc.cpp $(SRCdir)/YAML_Doc.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/YAML_Element.o: $(SRCdir)/YAML_Element.cpp $(SRCdir)/YAML_Element.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeDotProduct.o: $(SRCdir)/ComputeDotProduct.cpp $(SRCdir)/ComputeDotProduct.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeDotProduct_ref.o: $(SRCdir)/ComputeDotProduct_ref.cpp $(SRCdir)/ComputeDotProduct_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/finalize.o: $(SRCdir)/finalize.cpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/init.o: $(SRCdir)/init.cpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/mytimer.o: $(SRCdir)/mytimer.cpp $(SRCdir)/mytimer.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeOptimalShapeXYZ.o: $(SRCdir)/ComputeOptimalShapeXYZ.cpp $(SRCdir)/ComputeOptimalShapeXYZ.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeSPMV.o: $(SRCdir)/ComputeSPMV.cpp $(SRCdir)/ComputeSPMV.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeSPMV_ref.o: $(SRCdir)/ComputeSPMV_ref.cpp $(SRCdir)/ComputeSPMV_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeSYMGS.o: $(SRCdir)/ComputeSYMGS.cpp $(SRCdir)/ComputeSYMGS.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeSYMGS_ref.o: $(SRCdir)/ComputeSYMGS_ref.cpp $(SRCdir)/ComputeSYMGS_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeWAXPBY.o: $(SRCdir)/ComputeWAXPBY.cpp $(SRCdir)/ComputeWAXPBY.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeWAXPBY_ref.o: $(SRCdir)/ComputeWAXPBY_ref.cpp $(SRCdir)/ComputeWAXPBY_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeMG_ref.o: $(SRCdir)/ComputeMG_ref.cpp $(SRCdir)/ComputeMG_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeMG.o: $(SRCdir)/ComputeMG.cpp $(SRCdir)/ComputeMG.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeProlongation_ref.o: $(SRCdir)/ComputeProlongation_ref.cpp $(SRCdir)/ComputeProlongation_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/ComputeRestriction_ref.o: $(SRCdir)/ComputeRestriction_ref.cpp $(SRCdir)/ComputeRestriction_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/GenerateCoarseProblem.o: $(SRCdir)/GenerateCoarseProblem.cpp $(SRCdir)/GenerateCoarseProblem.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/CheckAspectRatio.o: $(SRCdir)/CheckAspectRatio.cpp $(SRCdir)/CheckAspectRatio.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

src/OutputFile.o: $(SRCdir)/OutputFile.cpp $(SRCdir)/OutputFile.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -I$(INCdir) $< -o $@

