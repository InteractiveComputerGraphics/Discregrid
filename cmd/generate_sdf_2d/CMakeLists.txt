# Eigen library.
find_package(Eigen3 REQUIRED)

# Set include directories.
include_directories(
	../../extern
	../../discregrid/include
	${EIGEN3_INCLUDE_DIR}
)

if(WIN32)
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_USE_MATH_DEFINES)
endif(WIN32)

# Enable PolygonDistance to use clipper2 reliant constructor.
if (USE_CLIPPER2)
	add_definitions(-DIS_CLIPPER_ENABLED)
endif()

# OpenMP support.
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
	if (CMAKE_VERSION VERSION_GREATER "3.8")
		link_libraries(OpenMP::OpenMP_CXX)
	else()
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
endif()

add_executable(GenerateSDF2D
	main.cpp
)

add_dependencies(GenerateSDF2D
	Discregrid
)

target_link_libraries(GenerateSDF2D
	Discregrid
	Clipper2
	Clipper2utils
)

set_target_properties(GenerateSDF2D PROPERTIES FOLDER Cmd)
