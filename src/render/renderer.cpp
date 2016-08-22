#include "renderer.h"
#include "scalar_range_scaler.h"

#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkExtractEdges.h>


Renderer::Renderer() :
    vtk_renderer( vtkSmartPointer<vtkRenderer>::New() ),
    reader( vtkSmartPointer<vtkXMLUnstructuredGridReader>::New() ),
    soln_surface_filter( vtkSmartPointer<vtkDataSetSurfaceFilter>::New() ),
    solution( new SolutionPipeline(soln_surface_filter) )
{
    reader->ReadFromInputStringOn();
    // Apply surface filter so can convert to polydata
    soln_surface_filter->SetInputConnection(reader->GetOutputPort());
}

void Renderer::load_solution(const std::string& vtu)
{
    reader->SetInputString(vtu);
    if (initialize) { // setup that can only happen after setting input string
        add_actors();
        initialize = false;
    }
}

vtkSmartPointer<vtkRenderer> Renderer::get_vtk_renderer()
{
    return vtk_renderer;
}

std::vector<std::string> Renderer::get_dataset_labels()
{
    return solution->get_dataset_labels();
}

void Renderer::set_displayed_dataset(const std::string &label)
{
    solution->set_displayed_dataset(label);
}

void Renderer::add_actors()
{
    vtk_renderer->AddActor(solution->get_actor());

    auto grid = create_grid_actor();
    vtk_renderer->AddActor(grid);
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


// todo https://github.com/OpenChemistry/avogadrolibs/blob/master/utilities/vtktesting/imageregressiontest.h
