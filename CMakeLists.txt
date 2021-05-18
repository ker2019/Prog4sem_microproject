cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(microproject)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
	IMPORTED_LOCATION "/home/sandbox/.local/lib/python3.8/site-packages/gmsh-git-Linux64-sdk/lib/libgmsh.so"
	INTERFACE_INCLUDE_DIRECTORIES "/home/sandbox/.local/lib/python3.8/site-packages/gmsh-git-Linux64-sdk/include"
)
set(GMSH_LIBRARY gmsh)

set(VTK_DIR "/home/sandbox/vtkbase/lib/cmake/vtk-7.1")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#add_executable(wave2 wave2.cpp wave2.h)
add_executable(envelope envelope.cpp envelope.h)

add_executable(meshgen meshgen.cpp geometry.cpp geometry.hpp)
target_link_libraries(meshgen dolfin ${GMSH_LIBRARY})

#target_link_libraries(wave2 dolfin)
target_link_libraries(envelope dolfin ${VTK_LIBRARIES})
