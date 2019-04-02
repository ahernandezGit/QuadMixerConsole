QT -= gui

CONFIG += c++11 console
CONFIG -= app_bundle

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS
QMAKE_CXXFLAGS += -std=c++0x
# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

EIGEN_PATH  =  /usr/include/eigen3
VCGLIB_PATH = $$PWD/vcglib/
LIBIGL_PATH = $$PWD/igl/
BUILDDIR = $$PWD/build/
PRECOMPILED_HEADER =$$PWD/precompiled_libraries.h
CONFIG +=precompiled_header
CONFIG += no_keywords

INCLUDEPATH += $$VCGLIB_PATH/ $$EIGEN_PATH $$LIBIGL_PATH
#LIBS += -L$$GUROBI_PATH/lib
#LIBS += -lgurobi_g++5.2 -lgurobi80
LIBS += -lboost_log -lboost_thread -lboost_iostreams -lboost_filesystem -lpthread -lboost_system -lCGAL -lgmp -lmpfr -fopenmp

# Directories
OBJECTS_DIR =   $$BUILDDIR/obj
MOC_DIR =       $$BUILDDIR/moc
RCC_DIR =       $$BUILDDIR/rccFind
UI_DIR =        $$BUILDDIR/ui
DESTDIR =       $$PWD/bin


SOURCES += \
        main.cpp \
    meshes.cpp \
    myutils.cpp \
    vectorcone.cpp \
    utils/smoothintersectioncurve.cpp \
    PatternsTakayama/ktmethod/patchgen/extradefinition.cpp \
    PatternsTakayama/ktmethod/patchgen/generate_topology.cpp \
    PatternsTakayama/ktmethod/patchgen/generate_subtopology.cpp \
    PatternsTakayama/ktmethod/patchgen/get_default_parameter.cpp \
    PatternsTakayama/ktmethod/patchgen/get_param_str.cpp \
    PatternsTakayama/ktmethod/patchgen/get_boundary_geometry.cpp \
    vcglib/wrap/ply/plylib.cpp

HEADERS += \
    meshes.h \
    PatternsTakayama/ktmethod/lp_solve/lp_lib.h \
    precompiled_libraries.h \
    myutils.h \
    vectorcone.h \
    meshtypes.h \
    patch3d.h \
    lsd.h \
    parameterizationlscm.h \
    cdt2d.h \
    quad_tracer.h \
    orient_faces.h \
    utils/smoothintersectioncurve.h \
    PatternsTakayama/patchg.h \
    PatternsTakayama/laplacianreconstruction.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_3_1.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_3_3.h \
    PatternsTakayama/ktmethod/patchgen/ilpgurobi.h \
    PatternsTakayama/ktmethod/patchgen/ILP.h \
    PatternsTakayama/ktmethod/patchgen/edgeloop.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_3_2.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_3_0.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_6_3.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_2_0.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_6_0.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_6_1.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_6_2.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_5_4.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_5_3.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_2_1.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_5_2.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_5_1.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_5_0.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_4_4.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_4_3.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_4_2.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_4_1.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_4_0.h \
    PatternsTakayama/ktmethod/patchgen/Pattern.h \
    PatternsTakayama/ktmethod/patchgen/generate_topology.h \
    PatternsTakayama/ktmethod/patchgen/Pattern_all.h \
    PatternsTakayama/ktmethod/patchgen/generate_subtopology.h \
    PatternsTakayama/ktmethod/patchgen/Permutation.h \
    PatternsTakayama/ktmethod/patchgen/decl.h \
    PatternsTakayama/ktmethod/patchgen/PatchParam.h \
    vcglib/wrap/ply/plylib.h \
    csvfile.h

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/PatternsTakayama/ktmethod/lp_solve/release/ -llpsolve55
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/PatternsTakayama/ktmethod/lp_solve/debug/ -llpsolve55
else:unix: LIBS += -L$$PWD/PatternsTakayama/ktmethod/lp_solve/ -llpsolve55

INCLUDEPATH += $$PWD/PatternsTakayama/ktmethod/lp_solve
DEPENDPATH += $$PWD/PatternsTakayama/ktmethod/lp_solve

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/PatternsTakayama/ktmethod/lp_solve/release/liblpsolve55.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/PatternsTakayama/ktmethod/lp_solve/debug/liblpsolve55.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/PatternsTakayama/ktmethod/lp_solve/release/lpsolve55.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/PatternsTakayama/ktmethod/lp_solve/debug/lpsolve55.lib
else:unix: PRE_TARGETDEPS += $$PWD/PatternsTakayama/ktmethod/lp_solve/liblpsolve55.a


