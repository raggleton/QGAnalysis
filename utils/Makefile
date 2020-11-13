LHAPDF=$(shell scram tool tag lhapdf LHAPDF_BASE)

plotPDFs: plotPDFs.cc
	g++ plotPDFs.cc -o plotPDFs.exe `${LHAPDF}/bin/lhapdf-config --cflags --ldflags` `root-config --cflags --libs`
