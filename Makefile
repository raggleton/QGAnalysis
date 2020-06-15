LIBRARY := SUHH2QGAnalysis
LHAPDFINC=$(shell scram tool tag lhapdf INCLUDE)
LHAPDFLIB=$(shell scram tool tag LHAPDF LIBDIR)

FASTJETINC=$(shell scram tool tag fastjet INCLUDE)
FASTJETLIB=$(shell scram tool tag fastjet LIBDIR)

FJCONTRIBINC=$(shell scram tool tag fastjet-contrib INCLUDE)
FJCONTRIBLIB=$(shell scram tool tag fastjet-contrib LIBDIR)
FJCONTRIBLIBNAME=$(shell scram tool tag fastjet-contrib LIB)

USERCXXFLAGS := -I${LHAPDFINC} -I${FASTJETINC} -I${FJCONTRIBINC} -Wno-unused-variable -Wall
# need to explicitly add libUnfold.so from ROOT, because they don't include it in
# their default from `root-config --libs`
# also see BuildFile.xml for Unfold lib
USERLDFLAGS := -lUnfold -lSUHH2core -lSUHH2common -lGenVector -lSUHH2JetMETObjects -L${LHAPDFLIB} -lLHAPDF -L${FASTJETLIB} -lfastjetplugins -lfastjettools -lfastjet -lsiscone -lsiscone_spherical -L${FJCONTRIBLIB} -l${FJCONTRIBLIBNAME}

# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
include ../Makefile.common
