Install cmake (version >= 2.8.8 required):
  wget http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz
  tar xf cmake-3.2.2.tar.gz
  cd cmake-3.2.2
  ./configure
  make
  sudo apt-get install checkinstall
  sudo checkinstall

Install Qt (from http://www.vtk.org/Wiki/VTK/Configure_and_Build#Qt_4.8..2A):
  mkdir qt5.5.1-install && cd qt5.5.1-install
  wget http://download.qt.io/official_releases/qt/5.5/5.5.1/qt-opensource-linux-x64-5.5.1.run
  chmod +x qt-opensource-linux-x64-5.5.1.run
  ./qt-opensource-linux-x64-5.5.1.run

Set up VTK (http://www.vtk.org/Wiki/VTK/Configure_and_Build#Download_VTK_Source_code):
  Install qt & vtk dependencies:
    sudo apt-get build-deps qt5-default
    sudo apt-get install build-essential libxt-dev libnetcdf-dev
  wget http://www.vtk.org/files/release/6.3/VTK-6.3.0.tar.gz
  tar xzf VTK-6.3.0.tar.gz
  mkdir VTK-6.3.0-build && cd VTK-6.3.0-build
  cmake -DVTK_QT_VERSION:STRING=5 \
-DQT_QMAKE_EXECUTABLE:PATH=/.../qt/qt5.5.1/5.5/gcc_64/bin/qmake \
-DVTK_Group_Qt:BOOL=ON \
-DCMAKE_PREFIX_PATH:PATH=/.../qt/qt5.5.1/5.5/gcc_64/lib/cmake/Qt5 \
-DBUILD_SHARED_LIBS:BOOL=ON ../VTK-6.3.0
 (note: you need to use the full path to DCMAKE_PREFIX_PATH, and possibly to QT_QMAKE_EXECUTABLE)
  make -j<# of cores>
  sudo checkinstall

Install deal.ii:
  tar -xvf deal.II-X.Y.Z.tar.gz
  mkdir dealii-build && cd dealii-build
  cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ../deal.II
  make install
  make test

Creating eclipse project:
  mkdir /home/liam/eclipse/workspace/fluid_flow/build
  cd /home/liam/eclipse/workspace/fluid_flow/build
  mkdir /home/liam/eclipse/workspace/fluid_flow/src
  cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_ECLIPSE_VERSION=4.5 ../src
Note: may need to upgrade cmake depending on the version of eclipse/cmake used 

Import project via:
  File->Import->General->Existing projects into workspace


For an overview of the vtk architecture see http://www.aosabook.org/en/vtk.html
