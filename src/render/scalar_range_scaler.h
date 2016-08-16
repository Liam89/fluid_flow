#ifndef SCALERRANGESCALER_H
#define SCALERRANGESCALER_H

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkMapper;

class ScalarRangeScaler : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(ScalarRangeScaler, vtkPolyDataAlgorithm)

  void PrintSelf(ostream& os, vtkIndent indent);

  static ScalarRangeScaler* New();

  void set_mapper(vtkSmartPointer<vtkMapper> new_mapper);

protected:
  ScalarRangeScaler();
  ~ScalarRangeScaler();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
  ScalarRangeScaler(const ScalarRangeScaler&);  // Not implemented.
  void operator=(const ScalarRangeScaler&);  // Not implemented.

  vtkSmartPointer<vtkMapper> mapper;
};

#endif // SCALERRANGESCALER_H
