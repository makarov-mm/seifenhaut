#include "MainWindow.h"
#include "Catenoid.h"
#include "FilmWidget.h"
#include <QCheckBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSlider>
#include <QVBoxLayout>

namespace {
constexpr double kSliderScale = 1000.0; // slider ticks per unit of h/R
}

MainWindow::MainWindow() {
    setWindowTitle("Seifenhaut — catenoid soap film");
    auto* central = new QWidget;
    auto* rootLayout = new QHBoxLayout(central);

    view_ = new FilmWidget;
    rootLayout->addWidget(view_, 1);

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
    sepSlider_->setRange(int(0.20 * kSliderScale), int(1.45 * kSliderScale));
    sepSlider_->setValue(int(0.80 * kSliderScale));
    sepValue_ = new QLabel("0.800");
    sepValue_->setFixedWidth(48);
    sepRow->addWidget(sepSlider_);
    sepRow->addWidget(sepValue_);
    p->addLayout(sepRow);

    auto* marks = new QLabel(
        QString("metastable above %1\nno equilibrium above %2")
            .arg(catenoid::rhoMetastable(), 0, 'f', 4)
            .arg(catenoid::rhoMax(), 0, 'f', 4));
    marks->setStyleSheet("color: gray;");
    p->addWidget(marks);
    p->addSpacing(12);

    stateLabel_ = new QLabel;
    neckLabel_ = new QLabel;
    areaLabel_ = new QLabel;
    for (QLabel* l : {stateLabel_, neckLabel_, areaLabel_}) {
        l->setTextFormat(Qt::RichText);
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

    rootLayout->addWidget(panel);
    setCentralWidget(central);

    connect(sepSlider_, &QSlider::valueChanged, this, [this](int v) {
        const double h = v / kSliderScale;
        sepValue_->setText(QString::number(h, 'f', 3));
        view_->setSeparation(h);
    });
    connect(dipButton_, &QPushButton::clicked, this,
            [this] { view_->reDip(); });
    connect(pauseBox, &QCheckBox::toggled, this,
            [this](bool on) { view_->setPaused(on); });
    connect(wireBox, &QCheckBox::toggled, this,
            [this](bool on) { view_->setWireframe(on); });
    connect(view_, &FilmWidget::statsChanged, this,
            &MainWindow::refreshStats);
}

void MainWindow::refreshStats() {
    const SoapFilm& f = view_->film();
    const double h = f.separation();
    const double R = f.ringRadius();
    const double rho = h / R;
    const bool exists = rho < catenoid::rhoMax();

    QString state, color;
    if (f.pinched()) {
        state = "collapsed — Goldschmidt solution (two disks)";
        color = "#e05050";
    } else if (!exists) {
        state = "collapsing — no equilibrium exists";
        color = "#e05050";
    } else if (rho > catenoid::rhoMetastable()) {
        state = "metastable — two disks already have less area";
        color = "#d0a030";
    } else {
        state = "stable — global area minimum";
        color = "#40b070";
    }
    stateLabel_->setText(
        QString("<b>State:</b> <span style='color:%1'>%2</span>")
            .arg(color, state));

    if (f.pinched()) {
        neckLabel_->setText("<b>Neck radius:</b> 0");
        areaLabel_->setText(
            QString("<b>Area:</b> %1 &nbsp;(2\u03C0R\u00B2 = %2)")
                .arg(f.totalArea(), 0, 'f', 4)
                .arg(2 * M_PI * R * R, 0, 'f', 4));
    } else {
        const QString analyticNeck =
            exists ? QString::number(catenoid::neckRadius(h, R), 'f', 4)
                   : QStringLiteral("\u2014");
        const QString analyticArea =
            exists ? QString::number(catenoid::area(h, R), 'f', 4)
                   : QStringLiteral("\u2014");
        neckLabel_->setText(QString("<b>Neck radius:</b> %1 &nbsp;(analytic %2)")
                                .arg(f.neckRadius(), 0, 'f', 4)
                                .arg(analyticNeck));
        areaLabel_->setText(QString("<b>Area:</b> %1 &nbsp;(analytic %2)")
                                .arg(f.totalArea(), 0, 'f', 4)
                                .arg(analyticArea));
    }
    dipButton_->setEnabled(f.pinched() && rho < catenoid::rhoMax());
}
