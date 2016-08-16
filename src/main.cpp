#include <QApplication>

#include "simulation/poisson.h"
#include "gui/mainwindow.h"

#include <string>
#include <iostream>
#include <fstream>

int main(int argc, char** argv)
{
    QApplication app{argc, argv};

    std::unique_ptr<Simulation::Poisson> simulation(new Simulation::Poisson());
    std::string ss = simulation->run();

    MainWindow mainWindow;
    mainWindow.set_result( ss );
    mainWindow.show();

    return app.exec();
}

