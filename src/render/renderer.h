#ifndef RENDERER_H
#define RENDERER_H

#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkAssignAttribute.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>

#include <string>
#include <vector>

class Renderer
{
public:
    Renderer();
    void load_solution(const std::string &vtu_contents);
    vtkSmartPointer<vtkRenderer> get_vtk_renderer();
    std::vector<std::string> get_data_labels();
private:
    vtkSmartPointer<vtkRenderer> vtk_renderer;
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader;
    vtkSmartPointer<vtkDataSetSurfaceFilter> soln_surface_filter;
    vtkSmartPointer<vtkAssignAttribute> soln_active_scalar;
    vtkSmartPointer<vtkDataSetMapper> soln_mapper;
    vtkSmartPointer<vtkActor> solution;
    std::vector<std::string> data_labels;

    void update_data_labels();
    void create_solution_actor();
    vtkSmartPointer<vtkActor> create_grid_actor();
    void add_actors();
};

#endif // RENDERER_H
