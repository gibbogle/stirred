VTKVER = 5

CONFIG      += uitools

contains(VTKVER,5) {
DEFINES += "VTK_VER=5"
INCLUDEPATH  += C:/VTK-VS10/include/vtk-5.10
}
contains(VTKVER,6) {
DEFINES += "VTK_VER=6"
INCLUDEPATH  += C:/VTK-6.3.0-VS10/include/vtk-6.3
}

message("VTKVER:")
message($$VTKVER)

INCLUDEPATH  += c:/qwt-5.2.1/src
INCLUDEPATH  += D:/ffmpeg/include

FORMS        += monolayer_GUI.ui plotwin.ui
HEADERS       = mainwindow.h qmylabel.h params.h plot.h log.h misc.h \
                result_set.h graphs.h field.h libmonolayer.h ImageSave.h \
                dialog.h myqgraphicsview.h qvideooutput.h qmycheckbox.h histogram_item.h \
                global.h drug.h plotwin.h qcustomplot.h
RESOURCES    += icons.qrc
SOURCES       = main.cpp mainwindow.cpp params.cpp plot.cpp misc.cpp \
                lognormal.cpp graphs.cpp field.cpp ImageSave.cpp \
                dialog.cpp myqgraphicsview.cpp  qvideooutput.cpp \
                information.cpp global.cpp histogram_item.cpp slots.cpp protocol.cpp \
                drug.cpp plotwin.cpp qcustomplot.cpp popup.cpp

# See cmake_link_command.txt for the full list of libraries that CMake links
contains(VTKVER,5) {
LIBS += -LC:/VTK-VS10/lib/vtk-5.10 -lQVTK -lvtkRendering -lvtkGraphics -lvtkImaging -lvtkIO -lvtkFiltering -lvtkCommon \
-lvtkpng -lvtktiff -lvtkjpeg -lvtkexpat -lvfw32 -lopengl32  \
-lwsock32 -lvtksys -lws2_32 -lvtkexoIIc -lvtkNetCDF \
-lvtklibxml2 -lvtkzlib -lvtkalglib \
-lgdi32 -lkernel32 -luser32 -lgdi32 -lwinspool -lshell32 -lole32 -loleaut32 -luuid -lcomdlg32 -ladvapi32
}
!contains(VTKVER,5) {
LIBS += -LC:/VTK-6.3.0-VS10/lib -lvtkGUISupportQT-6.3 -lvtkRenderingCore-6.3  -lvtkImagingCore-6.3 -lvtkIOCore-6.3 -lvtkCommonCore-6.3 \
-lvtkpng-6.3 -lvtktiff-6.3 -lvtkjpeg-6.3 -lvtkexpat-6.3 -lvtksys-6.3 -lvtkexoIIc-6.3 -lvtkNetCDF-6.3 -lvtklibxml2-6.3 -lvtkzlib-6.3 -lvtkalglib-6.3 \
-lvtkCommonColor-6.3 -lvtkCommonComputationalGeometry-6.3 \
-lvtkFiltersCore-6.3 -lvtkFiltersGeneral-6.3 -lvtkFiltersGeometry-6.3 -lvtkFiltersSources-6.3 -lvtkFiltersExtraction-6.3 \
-lvtkRenderingQT-6.3 -lvtkRenderingImage-6.3 -lvtkCommonExecutionModel-6.3 -lvtkCommonDataModel-6.3 -lvtkCommonMath-6.3 -lvtkCommonMisc-6.3 -lvtkCommonSystem-6.3 \
-lvtkCommonTransforms-6.3 -lvtkRenderingOpenGL-6.3\
-lvtkIOImage-6.3 -lvtkImagingGeneral-6.3 -lvtkImagingMath-6.3 -lvtkImagingColor-6.3 -lvtkImagingHybrid-6.3 -lvtkImagingMath-6.3 -lvtkImagingMorphological-6.3 -lvtkImagingSources-6.3 \
-lvtkInteractionImage-6.3 -lvtkInteractionStyle-6.3 -lvtkInteractionWidgets-6.3 \
-lws2_32 -lwsock32 -lvfw32 -lopengl32 -lgdi32 -lkernel32 -luser32 -lgdi32 -lwinspool -lshell32 -lole32 -loleaut32 -luuid -lcomdlg32 -ladvapi32
}

LIBS += -Lc:/qwt-5.2.1/lib -lqwt5
LIBS += C:/bin/monolayer.lib

QMAKE_LIBDIR += D:/ffmpeg/lib
win32:LIBS += avcodec.lib
win32:LIBS += avdevice.lib
win32:LIBS += avfilter.lib
win32:LIBS += avformat.lib
win32:LIBS += avutil.lib
win32:LIBS += postproc.lib
win32:LIBS += swresample.lib
win32:LIBS += swscale.lib

QT           += network

#QMAKE_CXXFLAGS += -Wno-deprecated -Wno-write-strings
#NOTE: fix for link with Qt-4.7
#QMAKE_LFLAGS    += -Wl,-enable-auto-import

# install
target.path = .
sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS icons 
sources.path = .
INSTALLS += target sources

DEFINES += __MSVS_DLL__
DEFINES += _CRT_SECURE_NO_WARNINGS
DEFINES += _CRT_SECURE_NO_DEPRECATE

QMAKE_LFLAGS += /OPT:NOREF
