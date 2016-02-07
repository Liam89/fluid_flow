#ifndef RENDERER_H
#define RENDERER_H

#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkActor.h>

#include <string>

class Renderer
{
public:
    Renderer(std::string vtu_contents);
    vtkSmartPointer<vtkRenderer> get_vtk_renderer();
private:
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader;
    vtkSmartPointer<vtkRenderer> vtk_renderer;

    vtkSmartPointer<vtkActor> create_solution_actor();
    vtkSmartPointer<vtkActor> create_grid_actor();
    void add_actors();
};

#endif // RENDERER_H
