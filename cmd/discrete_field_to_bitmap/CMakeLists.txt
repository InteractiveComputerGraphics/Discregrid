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
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(WIN32)

if ( CMAKE_COMPILER_IS_GNUCC )
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-multichar")
endif ( CMAKE_COMPILER_IS_GNUCC )

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

add_executable(DiscreteFieldToBitmap
	main.cpp
	bmp_file.hpp
	bmp_file.cpp
)

add_dependencies(DiscreteFieldToBitmap
	Discregrid
)

target_link_libraries(DiscreteFieldToBitmap
	Discregrid
)

set_target_properties(DiscreteFieldToBitmap PROPERTIES FOLDER Cmd)
