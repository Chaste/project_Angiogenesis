find_package(Chaste COMPONENTS cell_based)
chaste_do_project(Angiogenesis)

######### If building the Python modules we need to find boost python and python headers ###################### 
set(BUILD_ANGIOGENESIS_PYTHON "Build Python Bindings for Angiogenesis component" CACHE BOOL FALSE)
if(${BUILD_ANGIOGENESIS_PYTHON})
add_definitions(-DCHASTE_ANGIOGENESIS_PYTHON)
find_package(Boost COMPONENTS python REQUIRED)

# Uncomment for 'normal' way to find Python
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# If using virtualenv, need to manually set Python dirs due to find_package limitations
#include_directories("/home/grogan/.pyenv/versions/2.7.8/include/python2.7")
#set(PYTHON_LIBRARIES "/home/grogan/.pyenv/versions/2.7.8/lib/libpython2.7.so" )

######### Get the Chaste source files and libraries ###################### 
include_directories(${Chaste_INCLUDE_DIRS} ${Chaste_THIRD_PARTY_INCLUDE_DIRS})
set(PROJECT_ANGIO_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(PROJECT_ANGIO_LIB ${CMAKE_CURRENT_BINARY_DIR}/libchaste_project_Angiogenesis.so)
find_package(ProjectAngiogenesis PATHS ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${ProjectAngiogenesis_INCLUDE_DIRS})
endif(${BUILD_ANGIOGENESIS_PYTHON})

if(${BUILD_ANGIOGENESIS_PYTHON})

######### Build the Python modules ###################### 
set (ANGIOGENESIS_PYTHON_MODULES "")
set (ANGIOGENESIS_PYTHON_MODULE_LOCATIONS "")
list (APPEND ANGIOGENESIS_PYTHON_MODULES core)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/core/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES geometry)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/geometry/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES vessel)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/population/vessel/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES pde)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/pde/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES simulation)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/simulation/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES mesh)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/mesh/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES cell)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/population/cell/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES flow)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/simulation/flow/)
list (APPEND ANGIOGENESIS_PYTHON_MODULES angiogenesis)
list (APPEND ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/simulation/angiogenesis/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/src/python/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/)
list(LENGTH ANGIOGENESIS_PYTHON_MODULES len1)
math(EXPR len2 "${len1} - 1")
foreach(val RANGE ${len2})
    list(GET ANGIOGENESIS_PYTHON_MODULES ${val} python_module)
    list(GET ANGIOGENESIS_PYTHON_MODULE_LOCATIONS ${val} python_module_location)
    add_library(_${python_module} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/${python_module}.cpp)
    set_target_properties(_${python_module} PROPERTIES PREFIX "" LIBRARY_OUTPUT_DIRECTORY ${python_module_location} EXCLUDE_FROM_ALL FALSE)

    target_link_libraries(_${python_module} boost_python ${PYTHON_LIBRARIES} ${Chaste_THIRD_PARTY_LIBRARIES} ${Chaste_LIBRARIES} ${PROJECT_ANGIO_LIB})
endforeach()
add_custom_target(project_Angiogenesis_Python COMMAND make _core _geometry _vessel _pde _simulation _mesh _cell _flow _angiogenesis)
endif(${BUILD_ANGIOGENESIS_PYTHON})
