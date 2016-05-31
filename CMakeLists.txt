find_package(Chaste COMPONENTS cell_based)

######### If building the Python modules we need to find boost python and python headers ###################### 
find_package(Boost COMPONENTS python)

# Uncomment for 'normal' way to find Python
#find_package(PythonLibs REQUIRED)
#include_directories(${PYTHON_INCLUDE_DIRS})

# If using virtualenv, need to manually set Python dirs due to find_package limitations
include_directories("/home/grogan/.pyenv/versions/2.7.8/include/python2.7")
set(PYTHON_LIBRARIES "/home/grogan/.pyenv/versions/2.7.8/lib/libpython2.7.so" )

######### Get the Chaste source files and libraries ###################### 
include_directories(${Chaste_INCLUDE_DIRS} ${Chaste_THIRD_PARTY_INCLUDE_DIRS})
set(PROJECT_ANGIO_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(PROJECT_ANGIO_LIB ${CMAKE_CURRENT_BINARY_DIR}/libchaste_project_Angiogenesis.so)
find_package(ProjectAngiogenesis PATHS ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${ProjectAngiogenesis_INCLUDE_DIRS})

chaste_do_project(Angiogenesis)

######### Build the Python modules ###################### 
set (CASIE_TARGET_LIBS "")
list (APPEND CASIE_TARGET_LIBS core)
list (APPEND CASIE_TARGET_LIBS geometry)
list (APPEND CASIE_TARGET_LIBS vessel)
list (APPEND CASIE_TARGET_LIBS pde)
list (APPEND CASIE_TARGET_LIBS simulation)
list (APPEND CASIE_TARGET_LIBS mesh)
list (APPEND CASIE_TARGET_LIBS cell)
list (APPEND CASIE_TARGET_LIBS flow)
list (APPEND CASIE_TARGET_LIBS angiogenesis)

foreach(casie_lib ${CASIE_TARGET_LIBS})
    add_library(_${casie_lib} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/wrapper/${casie_lib}.cpp)
    set_target_properties(_${casie_lib} PROPERTIES PREFIX "" )
    target_link_libraries(_${casie_lib} boost_python ${PYTHON_LB_PATH} ${Chaste_THIRD_PARTY_LIBRARIES} ${Chaste_LIBRARIES} ${PROJECT_ANGIO_LIB})
endforeach()