#ifndef RENDERER_H
#define RENDERER_H

#include "solution_pipeline.h"

#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkActor.h>

#include <string>
#include <vector>
#include <memory>

class Renderer
{
public:
    Renderer();
    void load_solution(const std::string &vtu_contents);
    vtkSmartPointer<vtkRenderer> get_vtk_renderer();
    std::vector<std::string> get_dataset_labels();
    void set_displayed_dataset(const std::string &label);
private:
    vtkSmartPointer<vtkRenderer> vtk_renderer;
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader;
    vtkSmartPointer<vtkDataSetSurfaceFilter> soln_surface_filter;
    vtkSmartPointer<vtkAssignAttribute> soln_active_scalar;
    std::unique_ptr<SolutionPipeline> solution;
    bool initialize = true;

    vtkSmartPointer<vtkActor> create_grid_actor();
    void add_actors();
};

#endif // RENDERER_H
