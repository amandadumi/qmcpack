if (PKG_CONFIG_FOUND)
    message( STATUS "found pkg-config executable")

    # point pkg-config to the location of this tree's libint2.pc
    set(ENV{PKG_CONFIG_PATH} "${LAMMPS_ROOT}/:$ENV{PKG_CONFIG_PATH}")
    message( STATUS "$ENV{PKG_CONFIG_PATH}" )
    if(LAMMPS_FIND_QUIETLY)
        pkg_check_modules(PC_LAMMPS QUIET liblammps)
    else()
        pkg_check_modules(PC_LAMMPS liblammps)
    endif()
    set(LAMMPS_VERSION ${PC_LAMMPS_VERSION})
    message( STATUS "Lammps version is ${PC_LAMMPS_VERSION}" )

    find_path(LAMMPS_INCLUDE_DIR
            NAMES lammps.h
	    PATHS ${LAMMPS_ROOT}/../src/
            )
    message( STATUS "Lammps include dirs were ${PC_LAMMPS_INCLUDE_DIRS}" )
    message( STATUS "Lammps include dirs are ${LAMMPS_INCLUDE_DIR}" )

    # if (LAMMPS_SHARED_LIBRARY_ONLY)
    #     message( STATUS "Lammps SHARED LIBRARY ONLY" )
    #     set(_LAMMPS_LIB_NAMES "liblammps.so" "liblammps.dylib")
    # else (LAMMPS_SHARED_LIBRARY_ONLY)
    #     set(_LAMMPS_LIB_NAMES "liblammps")
    # endif(LAMMPS_SHARED_LIBRARY_ONLY)
    set(_LAMMPS_LIB_NAMES "liblammps.so" "liblammps.dylib")
    message( STATUS "Lammps iLIBRARY NAMES ${_LAMMPS_LIB_NAMES}" )

    message(STATUS "LAMMPS linker flags are ${PC_LAMMPS_LDFLAGS}")
    message(STATUS "LAMMPS cflags flags are ${PC_LAMMPS_CFLAGS}")
    #set(LAMMPS_CFLAGS ${PC_LAMMPS_CFLAGS})
    #set(LAMMPS_LDFLAGS ${PC_LAMMPS_LDFLAGS})

    find_library(LAMMPS_LIBRARY NAMES ${_LAMMPS_LIB_NAMES} HINTS ${LAMMPS_ROOT})
    message( STATUS "Lammps LIBRARY NAMES ${LAMMPS_LIBRARY}" )

    mark_as_advanced(LAMMPS_FOUND LAMMPS_INCLUDE_DIR LAMMPS_LIBRARY LAMMPS_VERSION LAMMPS_CFLAGS LAMMPS_LDFLAGS)
    mark_as_advanced(LAMMPS_FOUND LAMMPS_INCLUDE_DIR LAMMPS_LIBRARY LAMMPS_VERSION )

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(lammps
            FOUND_VAR LAMMPS_FOUND
            REQUIRED_VARS LAMMPS_INCLUDE_DIR
            VERSION_VAR LAMMPS_VERSION
            )

    if(LAMMPS_FOUND)
        set(LAMMPS_LIBRARIES ${LAMMPS_LIBRARY})
        set(LAMMPS_INCLUDE_DIRS ${LAMMPS_INCLUDE_DIR} ${LAMMPS_INCLUDE_DIR}/ML-SNAP)
	message( STATUS "Lammps include dirs on this linux machine ${LAMMPS_INCLUDE_DIR}" )
        # sanitize LAMMPS_INCLUDE_DIRS: remove duplicates and non-existent entries
        list(REMOVE_DUPLICATES LAMMPS_INCLUDE_DIRS)
        set(LAMMPS_INCLUDE_DIRS_SANITIZED )
        foreach(DIR IN LISTS LAMMPS_INCLUDE_DIRS)
            if (EXISTS ${DIR})
                list(APPEND LAMMPS_INCLUDE_DIRS_SANITIZED ${DIR})
            endif()
        endforeach()
        set(LAMMPS_INCLUDE_DIRS ${LAMMPS_INCLUDE_DIRS_SANITIZED})
    endif()

    if(LAMMPS_FOUND AND NOT TARGET Lammps::lammps)
    message( STATUS "we are here in creating the target!" )
    add_library(Lammps::lammps INTERFACE IMPORTED)
    set_target_properties(Lammps::lammps PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${LAMMPS_INCLUDE_DIRS}"
            )
    set_target_properties(Lammps::lammps PROPERTIES
            INTERFACE_LINK_LIBRARIES ${LAMMPS_LIBRARY}
            )
    set_target_properties(Lammps::lammps PROPERTIES
            INTERFACE_COMPILE_FEATURES "cxx_std_11"
            )
    endif()

else(PKG_CONFIG_FOUND)

    message(FATAL_ERROR "Could not find the required pkg-config executable")

endif(PKG_CONFIG_FOUND)
