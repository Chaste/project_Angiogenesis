# Build xerces

list(APPEND projects xerces)
set(_v 3)
set(xerces_version 3.1.${_v})
set(xerces_url "http://archive.apache.org/dist/xerces/c/3/sources/xerces-c-${xerces_version}.tar.gz")

set(xerces_cmds
CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
BUILD_COMMAND make all
INSTALL_COMMAND make install
)

ExternalProject_Add(xerces
  URL ${xerces_url}
  ${xerces_cmds}
)

ExternalProject_Get_Property(xerces install_dir)
set(XERCES_DIR "${install_dir}" CACHE INTERNAL "")
