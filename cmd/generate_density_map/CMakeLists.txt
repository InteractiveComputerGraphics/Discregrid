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

add_executable(GenerateDensityMap
	main.cpp
	gauss_quadrature.hpp
	gauss_quadrature.cpp
	sph_kernel.hpp
)

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

add_dependencies(GenerateDensityMap
	Discregrid
)

target_link_libraries(GenerateDensityMap
	Discregrid
)

set_target_properties(GenerateDensityMap PROPERTIES FOLDER Cmd)
