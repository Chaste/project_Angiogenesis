# Build OpenCascade

list(APPEND projects opencascade)

ExternalProject_Add(opencascade
  GIT_REPOSITORY https://github.com/tpaviot/oce.git
  GIT_TAG ""
  CMAKE_ARGS -DOCE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
  UPDATE_COMMAND ""
)

ExternalProject_Get_Property(opencascade install_dir)
set(OCE_DIR "${install_dir}" CACHE INTERNAL "")
