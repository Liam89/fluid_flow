CMake details:
  https://cmake.org/Wiki/Eclipse_CDT4_Generator
  https://www.dealii.org/8.3.0/index.html -> CMake in user projects


Creating eclipse project:
  mkdir /home/liam/eclipse/workspace/fluid_flow/build
  cd /home/liam/eclipse/workspace/fluid_flow/build
  mkdir /home/liam/eclipse/workspace/fluid_flow/src
  cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug ../src
Note: may need to upgrade cmake depending on the version of eclipse/cmake used 

Import project via:
  File->Import->General->Existing projects into workspace
