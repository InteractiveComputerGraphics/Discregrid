# - Try to find Discregrid library
#
# The find script can be invoked by using the following command:
#   find_package(Discregrid)
#
# Once done this will define
#
#  DISCREGRID_FOUND - System has found the discregrid library
#  DISCREGRID_INCLUDE_DIRS - Path to the discregrid include directory
#  DISCREGRID_LIBRARIES - Path to the static discregrid library

# ============================================================================
# _DISCREGRID_FIND_INCLUDE_DIR
# Internal function to find the include directories
#     _var = variable to set
#     _hdr = header file to look for
# ============================================================================
function(_DISCREGRID_FIND_INCLUDE_DIR _var _hdr)
    find_path(${_var} ${_hdr}
        PATHS
        $ENV{DISCREGRID_ROOT}
        ${DISCREGRID_ROOT}
	PATH_SUFFIXES
	    /include
    )

    if (${_var})
        set(DISCREGRID_INCLUDE_DIRS ${DISCREGRID_INCLUDE_DIRS} ${${_var}} PARENT_SCOPE)
        if (NOT DISCREGRID_SKIP_MARK_AS_ADVANCED)
            mark_as_advanced(${_var})
        endif()
    endif()
endfunction(_DISCREGRID_FIND_INCLUDE_DIR)

# ============================================================================
# _DISCREGRID_FIND_LIBRARY
# Internal function to find libraries packaged with DISCREGRID
#     _var = library variable to create
# ============================================================================
function(_DISCREGRID_FIND_LIBRARY _var _lib _mode)
    find_library(${_var}
	NAMES ${_lib}
	PATHS
        $ENV{DISCREGRID_ROOT}
        ${DISCREGRID_ROOT}
	PATH_SUFFIXES
	    /lib
    )

    if(${_var})
        set(DISCREGRID_LIBRARIES ${DISCREGRID_LIBRARIES} ${_mode} ${${_var}} PARENT_SCOPE)
        if(NOT DISCREGRID_SKIP_MARK_AS_ADVANCED)
            mark_as_advanced(${_var})
        endif()
    endif()
endfunction(_DISCREGRID_FIND_LIBRARY)

# ============================================================================
#
# main()
#
# ============================================================================

#
# Find all libraries and include directories.
#
_DISCREGRID_FIND_INCLUDE_DIR(DISCREGRID_DISCREGRIDH_INCLUDE_DIR Discregrid/All)
_DISCREGRID_FIND_LIBRARY(DISCREGRID_LIBRARY_RELEASE "discregrid" optimized)
_DISCREGRID_FIND_LIBRARY(DISCREGRID_LIBRARY_DEBUG "discregrid_d" debug)

#
# Try to enforce components.
#
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(DISCREGRID DEFAULT_MSG
    DISCREGRID_DISCREGRIDH_INCLUDE_DIR
    DISCREGRID_LIBRARY_DEBUG
    DISCREGRID_LIBRARY_RELEASE
)

if(NOT DISCREGRID_FOUND)
    set(DISCREGRID_INCLUDE_DIRS)
    set(DISCREGRID_LIBRARIES)
endif()

if(DISCREGRID_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES DISCREGRID_INCLUDE_DIRS)
endif()
