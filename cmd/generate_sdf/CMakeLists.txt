# Declare configuration file to embed supplied mesh files.
set(RESOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/resources/")
configure_file("../resource_path.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/resource_path.hpp")

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

add_executable(GenerateSDF
	main.cpp
	resource_path.hpp
)

add_dependencies(GenerateSDF
	Discregrid
)

target_link_libraries(GenerateSDF
	Discregrid
)

set_target_properties(GenerateSDF PROPERTIES FOLDER Cmd)
