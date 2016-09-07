#include <QApplication>

#include "simulation/poisson.h"
#include "simulation/ns_incompressible.h"
#include "gui/mainwindow.h"

#include <string>
#include <iostream>

int main(int argc, char** argv)
{
    QApplication app{argc, argv};

    // todo
    //std::unique_ptr<Simulation::Poisson> simulation(new Simulation::Poisson());
    std::unique_ptr<Simulation::SimulationBase<2>> simulation(new Simulation::NSIncompressible<2>());
    std::string ss = simulation->run();

    MainWindow mainWindow;
    mainWindow.set_result( ss );
    mainWindow.show();

    return app.exec();
}

