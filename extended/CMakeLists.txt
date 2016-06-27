# Copyright (c) 2005-2016, University of Oxford.
# All rights reserved.
# 
# University of Oxford means the Chancellor, Masters and Scholars of the
# University of Oxford, having an administrative office at Wellington
# Square, Oxford OX1 2JD, UK.
# 
# This file is part of Chaste.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  * Neither the name of the University of Oxford nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

cmake_minimum_required (VERSION 2.8.10)
project (Chaste_Angiogenesis_Extended)
set (Chaste_Angiogenesis_Extended_VERSION_MAJOR 3)
set (Chaste_Angiogenesis_Extended_VERSION_MINOR 5)

######### Build external projects ###################### 
set (ChasteAE_EXTERNAL_PROJECTS "")
list (APPEND ChasteAE_EXTERNAL_PROJECTS Boost)
list (APPEND ChasteAE_EXTERNAL_PROJECTS Xerces)
list (APPEND ChasteAE_EXTERNAL_PROJECTS Xsd)
list (APPEND ChasteAE_EXTERNAL_PROJECTS OpenMPI)
list (APPEND ChasteAE_EXTERNAL_PROJECTS PETSc)
list (APPEND ChasteAE_EXTERNAL_PROJECTS PETSc4Py)
list (APPEND ChasteAE_EXTERNAL_PROJECTS VTK)
list (APPEND ChasteAE_EXTERNAL_PROJECTS ITK)
list (APPEND ChasteAE_EXTERNAL_PROJECTS VMTK)
list (APPEND ChasteAE_EXTERNAL_PROJECTS FreeCad)
list (APPEND ChasteAE_EXTERNAL_PROJECTS Fenics)
list (APPEND ChasteAE_EXTERNAL_PROJECTS WxPython)
list (APPEND ChasteAE_EXTERNAL_PROJECTS Chaste)

foreach(ChasteAE_external ${ChasteAE_EXTERNAL_PROJECTS})
    option(BUILD_${ChasteAE_external} "Build ${ChasteAE_external}" ON)
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
