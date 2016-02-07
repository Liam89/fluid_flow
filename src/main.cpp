#include <QApplication>

#include "simulation/poisson.h"
#include "render/renderer.h"
#include "gui/mainwindow.h"

#include <string>
#include <iostream>

int main(int argc, char** argv)
{
    QApplication app{argc, argv};

    Simulation::Poisson poisson_problem;
    std::string ss = poisson_problem.run();
    Renderer renderer{ss};

    MainWindow mainWindow;
    mainWindow.set_renderer( renderer.get_vtk_renderer() );
    mainWindow.show();

    return app.exec();
}

