#pragma once
#include "SoapFilm.h"

#include <QMatrix4x4>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>

#include <vector>

// Renders the SoapFilm mesh and advances the relaxation within a fixed time
// budget per frame. The solver is deliberately owned by the widget: the mesh
// is small enough that stepping it on the GUI thread costs ~1 ms.
class FilmWidget : public QOpenGLWidget, protected QOpenGLFunctions_3_3_Core {
    Q_OBJECT
public:
    explicit FilmWidget(QWidget* parent = nullptr);
    ~FilmWidget() override;

    [[nodiscard]] const SoapFilm& film() const { return film_; }

    void setSeparation(double h);
    void reDip(); // fresh film at the current separation
    void setPaused(bool paused) { paused_ = paused; }
    void setWireframe(bool wireframe) { wireframe_ = wireframe; }

signals:
    void statsChanged(); // emitted once per frame after stepping

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void wheelEvent(QWheelEvent* e) override;

private:
    // Solver time budget per frame. 8 ms leaves half of a 60 Hz frame for
    // the buffer rebuild and rendering.
    static constexpr int kSimBudgetMs = 8;
    static constexpr int kMeshRows = 48, kMeshCols = 64;

    struct GpuVertex {
        float px, py, pz;
        float nx, ny, nz;
    };

    void advanceSimulation();
    void uploadFilm();
    void uploadRings();
    void upload(QOpenGLBuffer& vbo, QOpenGLVertexArrayObject& vao, const std::vector<GpuVertex>& vertices);

    SoapFilm film_;
    bool paused_ = false;
    bool wireframe_ = false;

    QOpenGLShaderProgram prog_;
    QOpenGLVertexArrayObject vaoFilm_, vaoRings_;
    QOpenGLBuffer vboFilm_{QOpenGLBuffer::VertexBuffer};
    QOpenGLBuffer vboRings_{QOpenGLBuffer::VertexBuffer};
    std::vector<GpuVertex> filmVertices_, ringVertices_;
    double ringsUploadedFor_ = -1.0; // separation the ring buffer was built for

    float yaw_ = 0.55f, pitch_ = 0.35f, dist_ = 3.2f;
    QPoint lastPos_;
    QMatrix4x4 proj_;
};
