# Build boost, based on https://github.com/OpenChemistry/openchemistry/blob/master/cmake/External_boost.cmake

list(APPEND projects boost)
set(_v 59)
set(boost_version 1.${_v}.0)
set(boost_url "http://sourceforge.net/projects/boost/files/boost/1.${_v}.0/boost_1_${_v}_0.tar.gz/download")

set(boost_with_args
  --with-serialization
  --with-filesystem
  --with-iostreams
  --with-program_options
  --with-system
  --with-thread
  --with-timer
  --with-regex
  --with-signals
  --with-python
)

set(boost_cmds
CONFIGURE_COMMAND ./bootstrap.sh --prefix=<INSTALL_DIR>
BUILD_COMMAND ./b2 ${boost_with_args}
INSTALL_COMMAND ./b2 ${boost_with_args} install
)

ExternalProject_Add(boost
  URL ${boost_url}
  ${boost_cmds}
  BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(boost install_dir)
set(BOOST_ROOT "${install_dir}" CACHE INTERNAL "")
