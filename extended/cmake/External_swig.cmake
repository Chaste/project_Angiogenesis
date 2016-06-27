# Build swig

list(APPEND projects swig)
set(swig_url "http://prdownloads.sourceforge.net/swig/swig-3.0.8.tar.gz")

set(swig_cmds
CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
BUILD_COMMAND make
INSTALL_COMMAND make install
)

ExternalProject_Add(swig
  URL ${swig_url}
  ${swig_cmds}
)

ExternalProject_Get_Property(swig install_dir)
set(SWIG_DIR "${install_dir}" CACHE INTERNAL "")

#ExternalProject_Add_Step(swig e1
#   COMMAND <SOURCE_DIR>/autogen.sh 
#   DEPENDERS configure
#   DEPENDEES patch
#)
