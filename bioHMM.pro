QT += core testlib
QT += gui

TARGET = bioHMM
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    pairhmm.cpp

HEADERS += \
    pairhmm.h \
    stl.h \
    math_utils.h \
    type.h \
    matrix.h \
    test/test_matrix.h \
    test/test_pairhmm.h \
    logsum.h \
    test/test_logsum.h

INCLUDEPATH += ./tools

LIBS += /Users/zhixingfeng/Documents/work/bioHMM/tools/UnitTest++/lib/libUnitTest++.dylib

