#pragma once
#include <QMainWindow>

class FilmWidget;
class QLabel;
class QPushButton;
class QSlider;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow();

private:
    QWidget* buildPanel();
    void refreshStats();

    FilmWidget* view_ = nullptr;
    QSlider* sepSlider_ = nullptr;
    QLabel* sepValue_ = nullptr;
    QLabel* stateLabel_ = nullptr;
    QLabel* neckLabel_ = nullptr;
    QLabel* areaLabel_ = nullptr;
    QPushButton* dipButton_ = nullptr;
};
