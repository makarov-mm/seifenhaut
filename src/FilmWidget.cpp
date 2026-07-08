#include "FilmWidget.h"
#include <QMouseEvent>
#include <QTimer>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {

const char* kVert = R"(#version 330 core
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
const char* kFrag = R"(#version 330 core
in vec3 vNormal;
in vec3 vWorldPos;
uniform vec3 uEye;
uniform vec4 uColor;
uniform int uUnlit;
out vec4 frag;
void main() {
    if (uUnlit == 1) { frag = uColor; return; }
    vec3 n = normalize(vNormal);
    vec3 v = normalize(uEye - vWorldPos);
    float ndv = abs(dot(n, v));
    vec3 l = normalize(vec3(0.4, 0.8, 0.5));
    float diff = 0.35 + 0.45 * abs(dot(n, l));
    float rim = pow(1.0 - ndv, 2.5);
    vec3 col = uColor.rgb * diff + vec3(0.85, 0.92, 1.0) * rim * 0.6;
    float alpha = uColor.a * (0.35 + 0.65 * rim);
    frag = vec4(col, alpha);
}
)";

} // namespace

FilmWidget::FilmWidget(QWidget* parent) : QOpenGLWidget(parent) {
    film_ = std::make_unique<SoapFilm>(48, 64, 1.0, 0.80);
    auto* timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, [this] {
        if (!paused_) {
            // Time-budgeted stepping: run as many solver steps as fit into
            // the budget, so the simulation runs at maximum speed on fast
            // machines without ever starving the render loop on slow ones.
            QElapsedTimer clock;
            clock.start();
            int done = 0;
            while (clock.nsecsElapsed() < simBudgetMs * 1'000'000LL &&
                   done < 5000) {
                film_->relax(20);
                done += 20;
            }
        }
        emit statsChanged();
        update();
    });
    timer->start(16);
    setMinimumSize(480, 420);
}

FilmWidget::~FilmWidget() {
    makeCurrent();
    vboFilm_.destroy();
    vboRings_.destroy();
    doneCurrent();
}

void FilmWidget::setSeparation(double h) { film_->setSeparation(h); }

void FilmWidget::reDip() {
    const double h = film_->separation();
    film_ = std::make_unique<SoapFilm>(48, 64, 1.0, h);
}

void FilmWidget::initializeGL() {
    initializeOpenGLFunctions();
    glClearColor(0.055f, 0.065f, 0.09f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kVert);
    prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kFrag);
    prog_.link();
    vaoFilm_.create();
    vaoRings_.create();
    vboFilm_.create();
    vboRings_.create();
}

void FilmWidget::resizeGL(int w, int h) {
    proj_.setToIdentity();
    proj_.perspective(38.0f, float(w) / float(std::max(1, h)), 0.05f, 50.0f);
}

void FilmWidget::rebuildFilmBuffer() {
    // Triangle soup with smooth (area-weighted accumulated) vertex normals,
    // rebuilt every frame — the mesh is a few thousand triangles, upload is
    // far cheaper than any dirty-tracking would be worth.
    std::vector<float> data;
    filmVertexCount_ = 0;
    for (const auto& p : film_->patches()) {
        std::vector<Vec3> nrm(p.X.size());
        for (const auto& t : p.tris) {
            const Vec3 n =
                (p.X[t[1]] - p.X[t[0]]).cross(p.X[t[2]] - p.X[t[0]]);
            nrm[t[0]] += n;
            nrm[t[1]] += n;
            nrm[t[2]] += n;
        }
        for (auto& n : nrm) {
            const double l = n.norm();
            if (l > 1e-14) n = n * (1.0 / l);
        }
        data.reserve(data.size() + p.tris.size() * 18);
        for (const auto& t : p.tris)
            for (int k = 0; k < 3; ++k) {
                const Vec3& v = p.X[t[k]];
                const Vec3& n = nrm[t[k]];
                data.insert(data.end(),
                            {float(v.x), float(v.y), float(v.z),
                             float(n.x), float(n.y), float(n.z)});
                ++filmVertexCount_;
            }
    }
    vaoFilm_.bind();
    vboFilm_.bind();
    vboFilm_.allocate(data.data(), int(data.size() * sizeof(float)));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 24, nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 24,
                          reinterpret_cast<void*>(12));
    vaoFilm_.release();
}

void FilmWidget::rebuildRingBuffer() {
    const int seg = 96;
    const double R = film_->ringRadius();
    const double h = film_->separation();
    std::vector<float> data;
    ringVertexCount_ = 0;
    for (double x : {-h / 2, h / 2})
        for (int j = 0; j < seg; ++j) {
            for (int k = 0; k < 2; ++k) {
                const double a = 2.0 * M_PI * (j + k) / seg;
                data.insert(data.end(),
                            {float(x), float(R * std::cos(a)),
                             float(R * std::sin(a)), 1.f, 0.f, 0.f});
                ++ringVertexCount_;
            }
        }
    vaoRings_.bind();
    vboRings_.bind();
    vboRings_.allocate(data.data(), int(data.size() * sizeof(float)));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 24, nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 24,
                          reinterpret_cast<void*>(12));
    vaoRings_.release();
}

void FilmWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    rebuildFilmBuffer();
    rebuildRingBuffer();

    QMatrix4x4 view;
    const QVector3D eye(dist_ * std::cos(pitch_) * std::sin(yaw_),
                        dist_ * std::sin(pitch_),
                        dist_ * std::cos(pitch_) * std::cos(yaw_));
    view.lookAt(eye, {0, 0, 0}, {0, 1, 0});
    QMatrix4x4 model; // world == model space; the film is centred already
    const QMatrix4x4 mvp = proj_ * view * model;

    prog_.bind();
    prog_.setUniformValue("uMvp", mvp);
    prog_.setUniformValue("uModel", model);
    prog_.setUniformValue("uEye", eye);

    // Wire rings: opaque, depth-written.
    prog_.setUniformValue("uUnlit", 1);
    prog_.setUniformValue("uColor", QVector4D(0.85f, 0.72f, 0.35f, 1.0f));
    vaoRings_.bind();
    glDrawArrays(GL_LINES, 0, ringVertexCount_);
    vaoRings_.release();

    // Film: translucent, two-sided, no depth writes so both layers show.
    prog_.setUniformValue("uUnlit", 0);
    prog_.setUniformValue("uColor", QVector4D(0.35f, 0.55f, 0.95f, 0.85f));
    glDepthMask(GL_FALSE);
    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, wireframe_ ? GL_LINE : GL_FILL);
    vaoFilm_.bind();
    glDrawArrays(GL_TRIANGLES, 0, filmVertexCount_);
    vaoFilm_.release();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDepthMask(GL_TRUE);
    prog_.release();
}

void FilmWidget::mousePressEvent(QMouseEvent* e) { lastPos_ = e->pos(); }

void FilmWidget::mouseMoveEvent(QMouseEvent* e) {
    const QPoint d = e->pos() - lastPos_;
    lastPos_ = e->pos();
    yaw_ += d.x() * 0.008f;
    pitch_ = std::clamp(pitch_ + d.y() * 0.008f, -1.45f, 1.45f);
    update();
}

void FilmWidget::wheelEvent(QWheelEvent* e) {
    dist_ = std::clamp(dist_ * std::pow(0.999f, float(e->angleDelta().y())),
                       1.2f, 12.0f);
    update();
}
