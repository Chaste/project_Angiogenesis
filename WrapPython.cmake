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

# Build the Python bindings
add_definitions(-DCHASTE_ANGIOGENESIS_PYTHON)
find_package(Boost COMPONENTS python REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${Chaste_INCLUDE_DIRS} ${Chaste_THIRD_PARTY_INCLUDE_DIRS})
set(PROJECT_ANGIO_LIB ${CMAKE_CURRENT_BINARY_DIR}/libchaste_project_Angiogenesis.so)

######### Build the Python modules ###################### 
set (ANGIOGENESIS_PYTHON_MODULES "")
set (ANGIOGENESIS_PYTHON_MODULE_LOCATIONS "")
list (APPEND ANGIOGENESIS_PYTHON_MODULES core)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/angiogenesis/core/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES geometry)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/angiogenesis/geometry/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES vessel)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/angiogenesis/population/vessel/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES pde)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/angiogenesis/pde/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES simulation)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/angiogenesis/simulation/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES mesh)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/angiogenesis/mesh/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES cell)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/angiogenesis/population/cell/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES flow)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/angiogenesis/simulation/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES angiogenesis)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/angiogenesis/simulation/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/src/python/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/ PATTERN "*.so" EXCLUDE)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test/python/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/test/)
list(LENGTH ANGIOGENESIS_PYTHON_MODULES len1)
math(EXPR len2 "${len1} - 1")
foreach(val RANGE ${len2})
    list(GET ANGIOGENESIS_PYTHON_MODULES ${val} python_module)
    list(GET ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${val} python_module_location)
    add_library(_chaste_project_Angiogenesis_${python_module} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/${python_module}.cpp)
    set_target_properties(_chaste_project_Angiogenesis_${python_module} PROPERTIES PREFIX "" LIBRARY_OUTPUT_DIRECTORY ${python_module_location})
    target_link_libraries(_chaste_project_Angiogenesis_${python_module} boost_python ${PYTHON_LIBRARIES} ${Chaste_THIRD_PARTY_LIBRARIES} ${Chaste_LIBRARIES} ${PROJECT_ANGIO_LIB})
endforeach()
add_custom_target(project_Angiogenesis_Python)
add_dependencies(project_Angiogenesis_Python 
    _chaste_project_Angiogenesis_core 
    _chaste_project_Angiogenesis_geometry 
    _chaste_project_Angiogenesis_vessel 
    _chaste_project_Angiogenesis_pde 
    _chaste_project_Angiogenesis_simulation 
    _chaste_project_Angiogenesis_mesh 
    _chaste_project_Angiogenesis_cell 
    _chaste_project_Angiogenesis_flow 
    _chaste_project_Angiogenesis_angiogenesis)
