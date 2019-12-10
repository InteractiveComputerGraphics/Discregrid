
set(HEADERS
	include/Discregrid/discrete_grid.hpp
	include/Discregrid/cubic_lagrange_discrete_grid.hpp
)

set(HEADERS_ACCELERATION

	include/Discregrid/acceleration/bounding_sphere_hierarchy.hpp
	include/Discregrid/acceleration/bounding_sphere.hpp
	include/Discregrid/acceleration/kd_tree.hpp
	include/Discregrid/acceleration/kd_tree.inl
)

set(HEADERS_MESH
	include/Discregrid/mesh/triangle_mesh.hpp
	include/Discregrid/mesh/entity_containers.hpp
	include/Discregrid/mesh/entity_iterators.hpp
	include/Discregrid/mesh/halfedge.hpp
)

set(HEADERS_GEOMETRY
	include/Discregrid/geometry/mesh_distance.hpp

	src/geometry/point_triangle_distance.hpp
)

set(HEADERS_UTILITY
	include/Discregrid/utility/serialize.hpp
	include/Discregrid/utility/lru_cache.hpp

	src/utility/timing.hpp
	src/utility/spinlock.hpp
)

set(SOURCES
	src/discrete_grid.cpp
	src/cubic_lagrange_discrete_grid.cpp
)

set(SOURCES_DATA
)

set(SOURCES_ACCELERATION
	src/acceleration/bounding_sphere_hierarchy.cpp
)

set(SOURCES_MESH
	src/mesh/entity_containers.cpp
	src/mesh/entity_iterators.cpp
	src/mesh/triangle_mesh.cpp
)

set(SOURCES_GEOMETRY
	src/geometry/mesh_distance.cpp
	src/geometry/point_triangle_distance.cpp
)

set(SOURCES_UTILITY
	src/utility/timing.cpp
)

macro(SOURCEGROUP name)
	string(TOLOWER ${name} name_lower)
	string(SUBSTRING ${name_lower} 0 1 FIRST_LETTER)
	string(TOUPPER ${FIRST_LETTER} FIRST_LETTER)
	string(REGEX REPLACE "^.(.*)" "${FIRST_LETTER}\\1" dname "${name_lower}")
	source_group("Source Files ${dname}" FILES ${SOURCES_${name}})
	source_group("Header Files ${dname}" FILES ${HEADERS_${name}})
endmacro()
SOURCEGROUP(MESH)
SOURCEGROUP(ACCELERATION)
SOURCEGROUP(DATA)
SOURCEGROUP(GEOMETRY)
SOURCEGROUP(UTILITY)

# OpenMP support.
if(APPLE)
	include(PatchOpenMPApple)
else()
	find_package(OpenMP REQUIRED)
endif()

if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Eigen library.
find_package(Eigen3 REQUIRED)

# Set include directories.
include_directories(
	include
	include/Discregrid
	${EIGEN3_INCLUDE_DIR}
	../extern
)

# Disable stupid MSVC compiler warnings.
if(WIN32)
	add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)
	add_definitions(-DNOMINMAX)
	add_definitions(-D_USE_MATH_DEFINES)
endif(WIN32)

set(DISCREGRID_SOURCE_FILES
	${HEADERS}
	${SOURCES}
	${HEADERS_ACCELERATION}
	${SOURCES_ACCELERATION}
	${SOURCES_DATA}
	${HEADERS_MESH}
	${SOURCES_MESH}
	${HEADERS_GEOMETRY}
	${SOURCES_GEOMETRY}
	${HEADERS_UTILITY}
	${SOURCES_UTILITY}
)

if(BUILD_AS_SHARED_LIBS)
	add_library(Discregrid SHARED ${DISCREGRID_SOURCE_FILES})
else()
	add_library(Discregrid ${DISCREGRID_SOURCE_FILES})
endif()

# Discregrid has the following dependencies.
#add_dependencies(Discregrid
#)

# Set link libraries.
target_link_libraries(Discregrid
)

install(TARGETS Discregrid
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib)

install(DIRECTORY include
DESTINATION .
PATTERN "*.hpp")
