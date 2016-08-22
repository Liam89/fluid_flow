#ifndef SOLUTION_PIPELINE_H
#define SOLUTION_PIPELINE_H

#include <vtkSmartPointer.h>

#include <string>
#include <vector>

class vtkAssignAttribute;
class vtkDataSetMapper;
class vtkDataSetSurfaceFilter;
class vtkActor;

class SolutionPipeline
{
public:
    SolutionPipeline(vtkSmartPointer<vtkDataSetSurfaceFilter> soln_surface_filter);
    std::vector<std::string> get_dataset_labels();
    void set_displayed_dataset(const std::string &label);
    void update();
    vtkSmartPointer<vtkActor> get_actor();
private:
    void update_dataset_labels();
    void create_solution_actor();
    vtkSmartPointer<vtkDataSetSurfaceFilter> soln_surface_filter;
    vtkSmartPointer<vtkAssignAttribute> soln_active_scalar;
    vtkSmartPointer<vtkDataSetMapper> soln_mapper;
    vtkSmartPointer<vtkActor> solution;
    std::vector<std::string> data_labels;
};

#endif // SOLUTION_PIPELINE_H
