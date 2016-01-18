#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkExtractEdges.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include <vtkGeometryFilter.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkDataSetSurfaceFilter.h>

#include <QApplication>

#include "simulation/laplace.h"

#include "gui/mainwindow.h"

#include <string>
#include <iostream>

vtkSmartPointer<vtkRenderer> render(std::string vtu)
{
	  //read all the data from the file
	  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
	    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	  reader->ReadFromInputStringOn();
	  reader->SetInputString(vtu);
      reader->Update();

      // Convert to polydata
      vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
        vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
      surfaceFilter->SetInputConnection(reader->GetOutputPort());
      surfaceFilter->Update();
      vtkSmartPointer<vtkPolyData> polydata = surfaceFilter->GetOutput();

      polydata->GetPointData()->SetActiveScalars("solution");

      // Find min and max z
      double bounds[2];
      polydata->GetScalarRange(bounds);
      double minz = bounds[0];
      double maxz = bounds[1];

      // grid
//      vtkSmartPointer<vtkExtractEdges> extractEdges =
//          vtkSmartPointer<vtkExtractEdges>::New();
//      extractEdges->SetInputConnection(reader->GetOutputPort());
//      extractEdges->Update();

	  //Create a mapper and actor
	  vtkSmartPointer<vtkDataSetMapper> mapper =
        vtkSmartPointer<vtkDataSetMapper>::New();
      mapper->SetInputConnection(surfaceFilter->GetOutputPort());
      mapper->SetScalarRange(minz, maxz);

	  vtkSmartPointer<vtkActor> actor =
	    vtkSmartPointer<vtkActor>::New();
      actor->SetMapper(mapper);

      //Create a renderer, and add the actor to the scene
	  vtkSmartPointer<vtkRenderer> renderer =
        vtkSmartPointer<vtkRenderer>::New();
      renderer->AddActor(actor);

      return renderer;
}

int main(int argc, char** argv)
{
    QApplication app(argc, argv);

    Simulation::Simulation laplace_problem;
    std::string ss = laplace_problem.run();
    vtkSmartPointer<vtkRenderer> renderer = render(ss);

    MainWindow mainWindow;
    mainWindow.set_renderer(renderer);
    mainWindow.show();

    return app.exec();
}

