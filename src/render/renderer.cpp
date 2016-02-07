#include "renderer.h"

#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataSetMapper.h>

#include <vtkDataSetSurfaceFilter.h>
#include <vtkWarpScalar.h>
#include <vtkExtractEdges.h>


Renderer::Renderer(std::string vtu)
{
    vtk_renderer = vtkSmartPointer<vtkRenderer>::New();

    reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->ReadFromInputStringOn();

    reader->SetInputString(vtu);
    add_actors();
}

vtkSmartPointer<vtkRenderer> Renderer::get_vtk_renderer()
{
    return vtk_renderer;
}

vtkSmartPointer<vtkActor> Renderer::create_solution_actor()
{
    // Apply surface filter so can set the active data points to the simulation results
    //  and scale the solution
    auto surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputConnection(reader->GetOutputPort());
    surfaceFilter->Update();

    reader->Update(); // Need to force update of pipeline before converting to
                      //  data object (polydata)
    vtkSmartPointer<vtkPolyData> polydata = surfaceFilter->GetOutput();
    // The simulation results are named 'solution' in the vtu xml
    polydata->GetPointData()->SetActiveScalars("solution");

    // Set the z axis to the solution
    auto warpScalar = vtkSmartPointer<vtkWarpScalar>::New();
    warpScalar->SetInputConnection(surfaceFilter->GetOutputPort());
    warpScalar->SetScaleFactor(1);
    warpScalar->UseNormalOn();
    warpScalar->SetNormal(0, 0, 1);

    //Create a mapper and scale the solution
    double bounds[2];
    polydata->GetScalarRange(bounds);
    auto minz = bounds[0], maxz = bounds[1];

    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputConnection(warpScalar->GetOutputPort());
    mapper->SetScalarRange(minz, maxz);

    auto solution = vtkSmartPointer<vtkActor>::New();
    solution->SetMapper(mapper);

    return solution;
}

vtkSmartPointer<vtkActor> Renderer::create_grid_actor()
{
    auto edges_filter = vtkSmartPointer<vtkExtractEdges>::New();
    edges_filter->SetInputConnection(reader->GetOutputPort());

    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputConnection(edges_filter->GetOutputPort());

    auto grid = vtkSmartPointer<vtkActor>::New();
    grid->SetMapper(mapper);

    return grid;
}

void Renderer::add_actors()
{
    auto solution = create_solution_actor();
    vtk_renderer->AddActor(solution);

    auto grid = create_grid_actor();
    vtk_renderer->AddActor(grid);
}
