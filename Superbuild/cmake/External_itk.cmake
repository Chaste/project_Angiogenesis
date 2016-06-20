# Build itk

list(APPEND projects itk)

set(itk_cmds
CMAKE_ARGS -DBUILD_SHARED_LIBS:BOOL=ON 
    -DBUILD_EXAMPLES:BOOL=off 
    -DITK_WRAP_PYTHON:BOOL=ON 
    -DITK_USE_FLAT_DIRECTORY_INSTALL:BOOL=ON 
    -DITK_USE_REVIEW:BOOL=ON # important
    -DModule_ITKReview:BOOL=ON # not sure if needed, playing it safe
    -DVTK_DIR:PATH=${VTK_DIR} 
    -DBUILD_TESTING:BOOL=OFF
    -DModule_ITKVtkGlue:BOOL=ON
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
)

ExternalProject_Add(itk
  DEPENDS vtk
  GIT_REPOSITORY http://itk.org/ITK.git
  GIT_TAG v4.9.1
  ${itk_cmds}
  UPDATE_COMMAND ""
    
)

ExternalProject_Get_Property(itk binary_dir)
set(ITK_DIR "${binary_dir}" CACHE INTERNAL "")
