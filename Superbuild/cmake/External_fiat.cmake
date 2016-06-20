# Build UFL

list(APPEND projects fiat)

ExternalProject_Add(fiat
  GIT_REPOSITORY https://bitbucket.org/fenics-project/fiat.git
  GIT_TAG ""
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND python <SOURCE_DIR>/setup.py install
  INSTALL_COMMAND ""
  BUILD_IN_SOURCE 1
)
