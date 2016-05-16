#-------------------------------------------------
#
# Project created by QtCreator 2016-04-18T23:03:34
#
#-------------------------------------------------

QT       += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = liza
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    qcustomplot.cpp \
    Solve.cpp

HEADERS  += mainwindow.h \
    qcustomplot.h \
    functions.h

FORMS    += mainwindow.ui
