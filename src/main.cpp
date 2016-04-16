#include <QApplication>

#include "simulation/poisson.h"
#include "render/renderer.h"
#include "gui/mainwindow.h"

#include <string>
#include <iostream>

int main(int argc, char** argv)
{
    QApplication app{argc, argv};

    std::unique_ptr<Simulation::Poisson> simulation(new Simulation::Poisson());
    std::string ss = simulation->run();
    Renderer renderer{ss};

    MainWindow mainWindow;
    mainWindow.set_renderer( renderer.get_vtk_renderer() );
    mainWindow.show();

    return app.exec();
}

