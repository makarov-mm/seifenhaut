#include "MainWindow.h"
#include "Catenoid.h"
#include "FilmWidget.h"

#include <QCheckBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSlider>
#include <QVBoxLayout>

#include <numbers>

namespace {

constexpr double kSliderScale = 1000.0; // slider ticks per unit of h/R
constexpr double kMinSeparation = 0.20, kMaxSeparation = 1.45, kInitialSeparation = 0.80;

QString fixed(double v, int decimals) { return QString::number(v, 'f', decimals); }

} // namespace

MainWindow::MainWindow() {
    setWindowTitle(QStringLiteral("Seifenhaut — catenoid soap film"));

    view_ = new FilmWidget;
    auto* central = new QWidget;
    auto* rootLayout = new QHBoxLayout(central);
    rootLayout->addWidget(view_, 1);
    rootLayout->addWidget(buildPanel());
    setCentralWidget(central);

    connect(sepSlider_, &QSlider::valueChanged, this, [this](int v) {
        const double h = v / kSliderScale;
        sepValue_->setText(fixed(h, 3));
        view_->setSeparation(h);
    });
    connect(dipButton_, &QPushButton::clicked, this, [this] { view_->reDip(); });
    connect(view_, &FilmWidget::statsChanged, this, &MainWindow::refreshStats);
}

QWidget* MainWindow::buildPanel() {
    auto* panel = new QWidget;
    panel->setFixedWidth(300);
    auto* p = new QVBoxLayout(panel);

    auto* heading = new QLabel(
        "Two rings, a soap film between them.\n"
        "Pull the rings apart and watch the\n"
        "catenoid lose its footing.");
    heading->setWordWrap(true);
    p->addWidget(heading);
    p->addSpacing(8);

    p->addWidget(new QLabel("Ring separation h / R"));
    auto* sepRow = new QHBoxLayout;
    sepSlider_ = new QSlider(Qt::Horizontal);
    sepSlider_->setRange(int(kMinSeparation * kSliderScale), int(kMaxSeparation * kSliderScale));
    sepSlider_->setValue(int(kInitialSeparation * kSliderScale));
    sepValue_ = new QLabel(fixed(kInitialSeparation, 3));
    sepValue_->setFixedWidth(48);
    sepRow->addWidget(sepSlider_);
    sepRow->addWidget(sepValue_);
    p->addLayout(sepRow);

    auto* marks = new QLabel(QStringLiteral("metastable above %1\nno equilibrium above %2")
                                 .arg(fixed(catenoid::rhoMetastable(), 4), fixed(catenoid::rhoMax(), 4)));
    marks->setStyleSheet("color: gray;");
    p->addWidget(marks);
    p->addSpacing(12);

    stateLabel_ = new QLabel;
    neckLabel_ = new QLabel;
    areaLabel_ = new QLabel;
    for (QLabel* l : {stateLabel_, neckLabel_, areaLabel_}) {
        l->setTextFormat(Qt::RichText);
        l->setWordWrap(true);
        p->addWidget(l);
    }
    p->addSpacing(12);

    dipButton_ = new QPushButton("Re-dip film");
    auto* pauseBox = new QCheckBox("Pause relaxation");
    auto* wireBox = new QCheckBox("Wireframe");
    p->addWidget(dipButton_);
    p->addWidget(pauseBox);
    p->addWidget(wireBox);
    p->addStretch(1);
    connect(pauseBox, &QCheckBox::toggled, view_, &FilmWidget::setPaused);
    connect(wireBox, &QCheckBox::toggled, view_, &FilmWidget::setWireframe);

    auto* footer = new QLabel(
        "The film is a triangle mesh relaxed by\n"
        "monotone descent on surface area\n"
        "(surface tension acting alone).\n"
        "Pinch-off and the Goldschmidt state\n"
        "emerge from the flow — nothing here\n"
        "is scripted.");
    footer->setStyleSheet("color: gray;");
    footer->setWordWrap(true);
    p->addWidget(footer);
    return panel;
}

void MainWindow::refreshStats() {
    const SoapFilm& f = view_->film();
    const double h = f.separation();
    const double R = f.ringRadius();
    const double rho = h / R;
    const bool exists = rho < catenoid::rhoMax();

    const auto [state, color] = [&]() -> std::pair<QString, QString> {
        if (f.pinched()) return {"collapsed — Goldschmidt solution (two disks)", "#e05050"};
        if (!exists) return {"collapsing — no equilibrium exists", "#e05050"};
        if (rho > catenoid::rhoMetastable()) return {"metastable — two disks already have less area", "#d0a030"};
        return {"stable — global area minimum", "#40b070"};
    }();
    stateLabel_->setText(QStringLiteral("<b>State:</b> <span style='color:%1'>%2</span>").arg(color, state));

    if (f.pinched()) {
        neckLabel_->setText(QStringLiteral("<b>Neck radius:</b> 0"));
        areaLabel_->setText(QStringLiteral("<b>Area:</b> %1 &nbsp;(2πR² = %2)")
                                .arg(fixed(f.totalArea(), 4), fixed(2 * std::numbers::pi * R * R, 4)));
    } else {
        const QString dash = QStringLiteral("—");
        neckLabel_->setText(QStringLiteral("<b>Neck radius:</b> %1 &nbsp;(analytic %2)")
                                .arg(fixed(f.neckRadius(), 4), exists ? fixed(catenoid::neckRadius(h, R), 4) : dash));
        areaLabel_->setText(QStringLiteral("<b>Area:</b> %1 &nbsp;(analytic %2)")
                                .arg(fixed(f.totalArea(), 4), exists ? fixed(catenoid::area(h, R), 4) : dash));
    }
    dipButton_->setEnabled(f.pinched() && exists);
}
