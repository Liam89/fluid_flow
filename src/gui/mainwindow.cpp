#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <vtkRenderWindow.h>
#include <QVTKWidget.h>
#include <QString>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->vtkWidget->GetRenderWindow()->AddRenderer(renderer.get_vtk_renderer());
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::set_result(const std::string &vtu_contents)
{
    renderer.load_solution(vtu_contents);
}

void MainWindow::on_startButton_clicked()
{

}

void MainWindow::on_stopButton_clicked()
{

}
