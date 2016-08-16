#include "scaler_range_scaler.h"

#include <vtkSetGet.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkMapper.h>

vtkStandardNewMacro(ScalarRangeScaler)

ScalarRangeScaler::ScalarRangeScaler()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ScalarRangeScaler::~ScalarRangeScaler()
{
}

void ScalarRangeScaler::set_mapper(vtkSmartPointer<vtkMapper> new_mapper)
{
    mapper = new_mapper;
}

int ScalarRangeScaler::RequestData(vtkInformation* vtkNotUsed(request),
                                             vtkInformationVector **inputVector,
                                             vtkInformationVector *outputVector)
{
  vtkPolyData* input = vtkPolyData::GetData(inputVector[0],0);
  vtkPolyData* output = vtkPolyData::GetData(outputVector,0);

  auto bounds = input->GetScalarRange();
  mapper->SetScalarRange(bounds);

  output->ShallowCopy(input);

  return 1;
}


//----------------------------------------------------------------------------
void ScalarRangeScaler::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
