#pragma once
#include "SoapFilm.h"
#include <QElapsedTimer>
#include <QMatrix4x4>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <memory>

// Renders the SoapFilm mesh and advances the relaxation a fixed number of
// substeps per frame. The solver is deliberately owned by the widget: the
// mesh is small enough that stepping it on the GUI thread costs ~1 ms.
class FilmWidget : public QOpenGLWidget, protected QOpenGLFunctions_3_3_Core {
    Q_OBJECT
public:
    explicit FilmWidget(QWidget* parent = nullptr);
    ~FilmWidget() override;

    SoapFilm& film() { return *film_; }

    void setSeparation(double h);
    void reDip();               // fresh film at the current separation
    void setPaused(bool p) { paused_ = p; }
    void setWireframe(bool w) { wireframe_ = w; }
    // Solver time budget per frame, milliseconds. 8 ms leaves half of a
    // 60 Hz frame for buffer rebuild and rendering.
    int simBudgetMs = 8;

signals:
    void statsChanged();        // emitted once per frame after stepping

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void wheelEvent(QWheelEvent* e) override;

private:
    void rebuildFilmBuffer();
    void rebuildRingBuffer();

    std::unique_ptr<SoapFilm> film_;
    bool paused_ = false;
    bool wireframe_ = false;

    QOpenGLShaderProgram prog_;
    QOpenGLVertexArrayObject vaoFilm_, vaoRings_;
    QOpenGLBuffer vboFilm_{QOpenGLBuffer::VertexBuffer};
    QOpenGLBuffer vboRings_{QOpenGLBuffer::VertexBuffer};
    int filmVertexCount_ = 0;
    int ringVertexCount_ = 0;

    float yaw_ = 0.55f, pitch_ = 0.35f, dist_ = 3.2f;
    QPoint lastPos_;
    QMatrix4x4 proj_;
};
