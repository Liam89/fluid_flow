#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <vtkRenderWindow.h>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::set_renderer(vtkSmartPointer<vtkRenderer> new_renderer)
{
    renderer = new_renderer;
    ui->vtkWidget->GetRenderWindow()->AddRenderer(renderer);
}

void MainWindow::on_startButton_clicked()
{

}

void MainWindow::on_stopButton_clicked()
{

}
