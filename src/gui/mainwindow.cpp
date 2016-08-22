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
    update_dataset_labels();
}

void MainWindow::update_dataset_labels()
{
    auto variableSelection = ui->variableSelection;
    variableSelection->clear();
    auto labels = renderer.get_dataset_labels();
    for (auto label : labels) {
        variableSelection->addItem( QString{label.c_str()} );
    }
    if ( !labels.empty() ) {
        renderer.set_displayed_dataset( labels[0] );
        variableSelection->setCurrentIndex(0);
    }
}

void MainWindow::on_startButton_clicked()
{

}

void MainWindow::on_stopButton_clicked()
{

}

void MainWindow::on_variableSelection_currentTextChanged(const QString &label)
{
    renderer.set_displayed_dataset( label.toStdString() );
    ui->vtkWidget->GetRenderWindow()->Render();
}
