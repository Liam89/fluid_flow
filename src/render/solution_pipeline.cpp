#include "solution_pipeline.h"
#include "scalar_range_scaler.h"

#include <vtkAssignAttribute.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPointData.h>
#include <vtkDataSetMapper.h>
#include <vtkWarpScalar.h>
#include <vtkActor.h>

SolutionPipeline::SolutionPipeline(vtkSmartPointer<vtkDataSetSurfaceFilter> soln_surface_filter) :
    soln_surface_filter(soln_surface_filter),
    soln_active_scalar( vtkSmartPointer<vtkAssignAttribute>::New() ),
    soln_mapper( vtkSmartPointer<vtkDataSetMapper>::New() ),
    solution( vtkSmartPointer<vtkActor>::New() )
{
    soln_active_scalar->SetInputConnection(soln_surface_filter->GetOutputPort());
    create_solution_actor();
}

void SolutionPipeline::update()
{
    update_dataset_labels();
}

void SolutionPipeline::update_dataset_labels() //todo update on get_data_labels
{
    // todo only update if modified
    soln_active_scalar->Update(); // Need to force update of pipeline before converting to
                                   //  data object
    data_labels.clear();
    auto point_data = soln_active_scalar->GetPolyDataOutput()->GetPointData();
    auto num_arrays = point_data->GetNumberOfArrays();
    for (int i = 0; i < num_arrays; ++i) {
        auto data_label = point_data->GetArrayName(i);
        assert(data_label != NULL);
        data_labels.push_back(data_label);
    }

//    if ( !data_labels.empty() ) {
//        soln_active_scalar->Assign(data_labels[0].c_str(),
//            vtkDataSetAttributes::SCALARS, vtkAssignAttribute::POINT_DATA);
//    }
}

void SolutionPipeline::create_solution_actor()
{
    // Set the z axis to the solution
    auto warp_scalar = vtkSmartPointer<vtkWarpScalar>::New();
    warp_scalar->SetInputConnection(soln_active_scalar->GetOutputPort());
    warp_scalar->SetScaleFactor(1);
    warp_scalar->UseNormalOn();
    warp_scalar->SetNormal(0, 0, 1);

    auto scaler = vtkSmartPointer<ScalarRangeScaler>::New();
    scaler->SetInputConnection(warp_scalar->GetOutputPort());
    scaler->set_mapper(soln_mapper.Get());

    soln_mapper->SetInputConnection(scaler->GetOutputPort());
    solution->SetMapper(soln_mapper);
}

vtkSmartPointer<vtkActor> SolutionPipeline::get_actor()
{
    return solution;
}

std::vector<std::string> SolutionPipeline::get_dataset_labels()
{
    update_dataset_labels();
    return data_labels;
}

void SolutionPipeline::set_displayed_dataset(const std::string &label)
{
    soln_active_scalar->Update();
    auto active_scalars = soln_active_scalar->GetPolyDataOutput()->GetPointData()->GetScalars(); // todo nicer way?
    if ( active_scalars == nullptr || active_scalars->GetName() != label ) {
        // SetActiveScalars doesn't update the solution actor, so vtkAssignAttribute is used
        soln_active_scalar->Assign(label.c_str(),
                vtkDataSetAttributes::SCALARS, vtkAssignAttribute::POINT_DATA);
    }
}
