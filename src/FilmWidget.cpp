#include "FilmWidget.h"

#include <QElapsedTimer>
#include <QMouseEvent>
#include <QTimer>
#include <QVector4D>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace {

constexpr const char* kVert = R"(#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
uniform mat4 uMvp;
uniform mat4 uModel;
out vec3 vNormal;
out vec3 vWorldPos;
void main() {
    gl_Position = uMvp * vec4(aPos, 1.0);
    vNormal = mat3(uModel) * aNormal;
    vWorldPos = (uModel * vec4(aPos, 1.0)).xyz;
}
)";

// Thin-film look: soft two-sided diffuse plus a Fresnel-ish rim so the
// silhouette glows the way a real soap film catches the light.
constexpr const char* kFrag = R"(#version 330 core
in vec3 vNormal;
in vec3 vWorldPos;
uniform vec3 uEye;
uniform vec4 uColor;
uniform bool uUnlit;
out vec4 frag;
void main() {
    if (uUnlit) { frag = uColor; return; }
    vec3 n = normalize(vNormal);
    vec3 v = normalize(uEye - vWorldPos);
    float ndv = min(abs(dot(n, v)), 1.0); // rounding can push it past 1 -> pow(-eps) = NaN
    vec3 l = normalize(vec3(0.4, 0.8, 0.5));
    float diff = 0.35 + 0.45 * abs(dot(n, l));
    float rim = pow(1.0 - ndv, 2.5);
    vec3 col = uColor.rgb * diff + vec3(0.85, 0.92, 1.0) * rim * 0.6;
    float alpha = uColor.a * (0.35 + 0.65 * rim);
    frag = vec4(col, alpha);
}
)";

constexpr QVector4D kRingColor(0.85f, 0.72f, 0.35f, 1.0f);
constexpr QVector4D kFilmColor(0.35f, 0.55f, 0.95f, 0.85f);

} // namespace

FilmWidget::FilmWidget(QWidget* parent)
    : QOpenGLWidget(parent), film_(kMeshRows, kMeshCols, 1.0, 0.80) {
    auto* timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, [this] {
        if (!paused_) advanceSimulation();
        emit statsChanged();
        update();
    });
    timer->start(16);
    setMinimumSize(480, 420);
}

FilmWidget::~FilmWidget() {
    makeCurrent(); // GL objects can only be released with their context current
    vboFilm_.destroy();
    vboRings_.destroy();
    vaoFilm_.destroy();
    vaoRings_.destroy();
    doneCurrent();
}

void FilmWidget::setSeparation(double h) { film_.setSeparation(h); }

void FilmWidget::reDip() { film_ = SoapFilm(kMeshRows, kMeshCols, 1.0, film_.separation()); }

// Time-budgeted stepping: run as many solver steps as fit into the budget,
// so the simulation runs at maximum speed on fast machines without ever
// starving the render loop on slow ones.
void FilmWidget::advanceSimulation() {
    QElapsedTimer clock;
    clock.start();
    for (int done = 0; done < 5000 && clock.nsecsElapsed() < kSimBudgetMs * 1'000'000LL; done += 20)
        film_.relax(20);
}

void FilmWidget::initializeGL() {
    initializeOpenGLFunctions();
    glClearColor(0.055f, 0.065f, 0.09f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    // colour blends with the film's alpha; the framebuffer alpha stays 1, so
    // the widget never turns see-through when the window system composites it
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ZERO, GL_ONE);
    if (!prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kVert) ||
        !prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kFrag) || !prog_.link())
        qFatal("Shader build failed: %s", qPrintable(prog_.log()));
    vaoFilm_.create();
    vaoRings_.create();
    vboFilm_.create();
    vboFilm_.setUsagePattern(QOpenGLBuffer::StreamDraw);
    vboRings_.create();
}

void FilmWidget::resizeGL(int w, int h) {
    proj_.setToIdentity();
    proj_.perspective(38.0f, float(w) / float(std::max(1, h)), 0.05f, 50.0f);
}

// (Re)fills the buffer and points the VAO at it: position, then normal.
void FilmWidget::upload(QOpenGLBuffer& vbo, QOpenGLVertexArrayObject& vao, const std::vector<GpuVertex>& vertices) {
    static_assert(sizeof(GpuVertex) == 6 * sizeof(float));
    vao.bind();
    vbo.bind();
    vbo.allocate(vertices.data(), int(vertices.size() * sizeof(GpuVertex)));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GpuVertex), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GpuVertex),
                          reinterpret_cast<const void*>(offsetof(GpuVertex, nx)));
    vao.release();
}

void FilmWidget::uploadFilm() {
    // Triangle soup with smooth (area-weighted accumulated) vertex normals,
    // rebuilt every frame — the mesh is a few thousand triangles, upload is
    // far cheaper than any dirty-tracking would be worth.
    filmVertices_.clear();
    std::vector<Vec3> normals;
    for (const FilmPatch& p : film_.patches()) {
        normals.assign(p.X.size(), Vec3{});
        for (const auto& [a, b, c] : p.tris) {
            const Vec3 n = (p.X[b] - p.X[a]).cross(p.X[c] - p.X[a]);
            normals[a] += n;
            normals[b] += n;
            normals[c] += n;
        }
        for (Vec3& n : normals) {
            const double l = n.norm();
            if (l > 1e-14) n *= 1.0 / l;
        }
        filmVertices_.reserve(filmVertices_.size() + p.tris.size() * 3);
        for (const auto& tri : p.tris)
            for (int k : tri) {
                const Vec3& v = p.X[k];
                const Vec3& n = normals[k];
                filmVertices_.push_back({float(v.x), float(v.y), float(v.z), float(n.x), float(n.y), float(n.z)});
            }
    }
    upload(vboFilm_, vaoFilm_, filmVertices_);
}

void FilmWidget::uploadRings() {
    const double h = film_.separation();
    if (h == ringsUploadedFor_) return; // the rings only move with the slider
    ringsUploadedFor_ = h;

    constexpr int kSegments = 96;
    const double R = film_.ringRadius();
    ringVertices_.clear();
    ringVertices_.reserve(2 * kSegments * 2);
    for (double x : {-h / 2, h / 2})
        for (int j = 0; j < kSegments; ++j)
            for (int k : {0, 1}) { // line list: both ends of every segment
                const double a = 2.0 * std::numbers::pi * (j + k) / kSegments;
                ringVertices_.push_back({float(x), float(R * std::cos(a)), float(R * std::sin(a)), 1.f, 0.f, 0.f});
            }
    upload(vboRings_, vaoRings_, ringVertices_);
}

void FilmWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    uploadFilm();
    uploadRings();

    const QVector3D eye(dist_ * std::cos(pitch_) * std::sin(yaw_),
                        dist_ * std::sin(pitch_),
                        dist_ * std::cos(pitch_) * std::cos(yaw_));
    QMatrix4x4 view;
    view.lookAt(eye, {0, 0, 0}, {0, 1, 0});
    const QMatrix4x4 model; // world == model space; the film is centred already

    prog_.bind();
    prog_.setUniformValue("uMvp", proj_ * view * model);
    prog_.setUniformValue("uModel", model);
    prog_.setUniformValue("uEye", eye);

    // Wire rings: opaque, depth-written.
    prog_.setUniformValue("uUnlit", GLint(1));
    prog_.setUniformValue("uColor", kRingColor);
    vaoRings_.bind();
    glDrawArrays(GL_LINES, 0, GLsizei(ringVertices_.size()));
    vaoRings_.release();

    // Film: translucent, two-sided, no depth writes so both layers show.
    prog_.setUniformValue("uUnlit", GLint(0));
    prog_.setUniformValue("uColor", kFilmColor);
    glDepthMask(GL_FALSE);
    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, wireframe_ ? GL_LINE : GL_FILL);
    vaoFilm_.bind();
    glDrawArrays(GL_TRIANGLES, 0, GLsizei(filmVertices_.size()));
    vaoFilm_.release();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDepthMask(GL_TRUE);
    prog_.release();
}

void FilmWidget::mousePressEvent(QMouseEvent* e) { lastPos_ = e->pos(); }

void FilmWidget::mouseMoveEvent(QMouseEvent* e) {
    if (!(e->buttons() & Qt::LeftButton)) return;
    const QPoint d = e->pos() - lastPos_;
    lastPos_ = e->pos();
    yaw_ += float(d.x()) * 0.008f;
    pitch_ = std::clamp(pitch_ + float(d.y()) * 0.008f, -1.45f, 1.45f);
    update();
}

void FilmWidget::wheelEvent(QWheelEvent* e) {
    dist_ = std::clamp(dist_ * std::pow(0.999f, float(e->angleDelta().y())), 1.2f, 12.0f);
    update();
}
