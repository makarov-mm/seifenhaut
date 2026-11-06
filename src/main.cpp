#include "MainWindow.h"

#include <QApplication>
#include <QSurfaceFormat>

int main(int argc, char** argv) {
    // The shaders are GLSL 3.30 core. macOS only hands out a modern context
    // when a core profile is requested up front, so ask before any window
    // (and any QOpenGLWidget) exists.
    QSurfaceFormat format;
    format.setVersion(3, 3);
    format.setProfile(QSurfaceFormat::CoreProfile);
    format.setDepthBufferSize(24);
    format.setSamples(4); // the wire rings and the wireframe are thin lines
    QSurfaceFormat::setDefaultFormat(format);

    QApplication app(argc, argv);
    QApplication::setApplicationName(QStringLiteral("Seifenhaut"));
    MainWindow window;
    window.resize(1100, 680);
    window.show();
    return QApplication::exec();
}
