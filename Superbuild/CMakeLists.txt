cmake_minimum_required (VERSION 2.6)
project (ChasteVT)
set (ChasteVT_VERSION_MAJOR 3)
set (ChasteVT_VERSION_MINOR 5)

######### Build external projects ###################### 
set (ChasteVT_EXTERNAL_PROJECTS "")
list (APPEND ChasteVT_EXTERNAL_PROJECTS Boost)
list (APPEND ChasteVT_EXTERNAL_PROJECTS Xerces)
list (APPEND ChasteVT_EXTERNAL_PROJECTS Xsd)
list (APPEND ChasteVT_EXTERNAL_PROJECTS OpenMPI)
list (APPEND ChasteVT_EXTERNAL_PROJECTS PETSc)
list (APPEND ChasteVT_EXTERNAL_PROJECTS PETSc4Py)
list (APPEND ChasteVT_EXTERNAL_PROJECTS VTK)
list (APPEND ChasteVT_EXTERNAL_PROJECTS ITK)
list (APPEND ChasteVT_EXTERNAL_PROJECTS VMTK)
list (APPEND ChasteVT_EXTERNAL_PROJECTS FreeCad)
list (APPEND ChasteVT_EXTERNAL_PROJECTS Fenics)
list (APPEND ChasteVT_EXTERNAL_PROJECTS WxPython)
list (APPEND ChasteVT_EXTERNAL_PROJECTS Chaste)

foreach(ChasteVT_external ${ChasteVT_EXTERNAL_PROJECTS})
    option(BUILD_${ChasteVT_external} "Build ${ChasteVT_external}" ON)
endforeach()

# Set up for building external projects
include(ExternalProject)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set_directory_properties(PROPERTIES EP_BASE "${CMAKE_BINARY_DIR}/third_party")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/third_party/lib")

# Find Python
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

if(BUILD_Xerces)
    include(External_xerces)
endif()

if(BUILD_Xsd)
    include(External_xsd)
endif()

if(BUILD_OpenMPI)
    include(External_openmpi)
endif()

if(BUILD_Boost)
    include(External_boost)
endif()

if(BUILD_PETSc)
    if(CMAKE_BUILD_TYPE MATCHES RELEASE)
        include(External_petsc_opt)
    else()
        include(External_petsc)
    endif()

#    if(BUILD_PETSc4Py)
#        include(External_petsc4py)
#    endif()
endif()

if(BUILD_VTK)
    include(External_vtk)
endif()

if(BUILD_ITK)
    include(External_itk)
endif()

if(BUILD_VMTK)
    include(External_vmtk)
endif()

if(BUILD_FreeCad)
    include(External_swig)
    include(External_opencascade)
    include(External_freecad)
endif()

OPTION(BUILD_Fenics "Off for now" OFF)
if(BUILD_Fenics)
    if (NOT BUILD_FreeCad)
        include(External_swig)
    endif()
    include(External_ufl)
    include(External_instant)
    include(External_fiat)
    include(External_ffc)
    include(External_dolfin)
endif()

if(BUILD_Chaste)
    include(External_chaste)
endif()
