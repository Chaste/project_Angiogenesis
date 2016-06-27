# Build the project as a Chaste component with Chaste only dependencies.
find_package(Chaste COMPONENTS cell_based)
chaste_do_project(Angiogenesis)

set(BUILD_ANGIOGENESIS_PYTHON "Build Python Bindings for Angiogenesis component" CACHE BOOL TRUE)
if(${BUILD_ANGIOGENESIS_PYTHON})
    include(${CMAKE_CURRENT_SOURCE_DIR}/ProjectAngiogenesisPython.cmake)
endif${BUILD_ANGIOGENESIS_PYTHON})

