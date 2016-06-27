# Build vmtk

list(APPEND projects vmtk)

set(vmtk_cmds
CMAKE_ARGS -DBUILD_SHARED_LIBS:BOOL=ON 
    -DUSE_SYSTEM_ITK:BOOL=ON 
    -DUSE_SYSTEM_VTK:BOOL=ON 
    -DITK_DIR:PATH=${ITK_DIR} 
    -DVTK_DIR:PATH=${VTK_DIR}
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
)

ExternalProject_Add(vmtk
  DEPENDS itk
  GIT_REPOSITORY http://github.com/vmtk/vmtk.git
  GIT_TAG v1.2
  ${vmtk_cmds}
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""  
)

ExternalProject_Get_Property(vmtk binary_dir)
set(VMTK_DIR "${binary_dir}" CACHE INTERNAL "")
