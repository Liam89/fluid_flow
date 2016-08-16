#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <vtkSmartPointer.h>
#include <render/renderer.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void set_result(const std::string &vtu_contents);

private Q_SLOTS:
    void on_startButton_clicked();

    void on_stopButton_clicked();
    void on_variableSelection_currentTextChanged(const QString &label);

private:
    Ui::MainWindow *ui;
    Renderer renderer;

    void update_data_labels();
};

#endif // MAINWINDOW_H
