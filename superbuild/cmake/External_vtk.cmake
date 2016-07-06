# Build vtk

list(APPEND projects vtk)

set(vtk_cmds
CMAKE_ARGS -DBUILD_SHARED_LIBS:BOOL=ON 
    -DBUILD_TESTING:BOOL=OFF 
    -DCMAKE_C_FLAGS=-DGLX_GLXEXT_LEGACY # ubuntu 16 bug
    -DCMAKE_CXX_FLAGS=-DGLX_GLXEXT_LEGACY # ubuntu 16 bug
    -DVTK_WRAP_PYTHON:BOOL=ON 
    -DVTK_USE_SYSTEM_HDF5:BOOL=OFF 
    -DHDF5_DIR=${HDF5_DIR} 
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
)

ExternalProject_Add(vtk
  DEPENDS petsc
  GIT_REPOSITORY http://vtk.org/VTK.git
  GIT_TAG v5.10.1
  ${vtk_cmds}
  UPDATE_COMMAND ""
    
)

ExternalProject_Get_Property(vtk binary_dir)
set(VTK_DIR "${binary_dir}" CACHE INTERNAL "")
