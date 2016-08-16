#include "renderer.h"

#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkAssignAttribute.h>
#include <vtkDataSetMapper.h>
#include <vtkWarpScalar.h>
#include <vtkExtractEdges.h>


Renderer::Renderer() :
    vtk_renderer( vtkSmartPointer<vtkRenderer>::New() ),
    reader( vtkSmartPointer<vtkXMLUnstructuredGridReader>::New() ),
    soln_surface_filter( vtkSmartPointer<vtkDataSetSurfaceFilter>::New() ),
    soln_active_scalar( vtkSmartPointer<vtkAssignAttribute>::New() ),
    soln_mapper( vtkSmartPointer<vtkDataSetMapper>::New() ),
    solution( vtkSmartPointer<vtkActor>::New() )
{
    reader->ReadFromInputStringOn();
}

void Renderer::load_solution(const std::string& vtu)
{
    reader->SetInputString(vtu);

    // Convert to polydata
    // Apply surface filter so can set the active data points to the simulation results
    //  and scale the solution
    soln_surface_filter->SetInputConnection(reader->GetOutputPort());
    soln_active_scalar->SetInputConnection(soln_surface_filter->GetOutputPort());

    update_data_labels();

    add_actors();
}

vtkSmartPointer<vtkRenderer> Renderer::get_vtk_renderer()
{
    return vtk_renderer;
}

std::vector<std::string> Renderer::get_data_labels()
{
    return data_labels;
}

void Renderer::update_data_labels()
{
    soln_surface_filter->Update(); // Need to force update of pipeline before converting to
                                   //  data object
    auto point_data = soln_surface_filter->GetOutput()->GetPointData();
    auto num_arrays = point_data->GetNumberOfArrays();
    assert(num_arrays != 0);
    for (int i = 0; i < num_arrays; ++i) {
        auto data_label = point_data->GetArrayName(i);
        assert(data_label != NULL);
        data_labels.push_back(data_label);
    }

    soln_active_scalar->Assign(data_labels[0].c_str(),
            vtkDataSetAttributes::SCALARS, vtkAssignAttribute::POINT_DATA);
}


void Renderer::add_actors()
{
    create_solution_actor();
    vtk_renderer->AddActor(solution);

    auto grid = create_grid_actor();
    vtk_renderer->AddActor(grid);
}

void Renderer::create_solution_actor()
{
    // Set the z axis to the solution
    auto warp_scalar = vtkSmartPointer<vtkWarpScalar>::New();
    warp_scalar->SetInputConnection(soln_active_scalar->GetOutputPort());
    warp_scalar->SetScaleFactor(1);
    warp_scalar->UseNormalOn();
    warp_scalar->SetNormal(0, 0, 1);

    soln_mapper->SetInputConnection(warp_scalar->GetOutputPort());
    solution->SetMapper(soln_mapper);
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
