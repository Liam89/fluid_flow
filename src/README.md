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

Install gmock & gtest:
  mkdir gmock && cd gmock
  wget https://github.com/google/googlemock/archive/release-1.7.0.tar.gz
  tar xzf release-1.7.0.tar.gz
  mv googlemock-release-1.7.0 gmock-1.7.0
  rm release-1.7.0.tar.gz
  wget https://github.com/google/googletest/archive/release-1.7.0.tar.gz
  tar xzf release-1.7.0.tar.gz
  mv googletest-release-1.7.0/ gmock-1.7.0/gtest
  rm release-1.7.0.tar.gz




Useful links:
Cmake tutorial https://www.johnlamp.net/cmake-tutorial-1-getting-started.html
Overview of the vtk architecture http://www.aosabook.org/en/vtk.html
