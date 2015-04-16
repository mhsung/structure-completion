cmake_minimum_required (VERSION 2.6)
project (QtViewer)

# add our macro directory to cmake search path
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set (LIBRARY_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../lib CACHE PATH "The directory where the all library files can be found.")
set (OPENMESH_DIR ${LIBRARY_ROOT_PATH}/OpenMesh CACHE PATH "The directory where the OpenMesh files can be found.")

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${OPENMESH_DIR}/cmake)
set (CMAKE_DEBUG_POSTFIX "d")

# include our cmake files
include (ACGCommon)
include (ACGOutput)

if (WIN32)
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
      -D_CRT_SECURE_NO_WARNINGS
  )
endif ()


include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11
support. Please use a different C++ compiler.")
endif()


# ========================================================================
# Add bundle targets here
# ========================================================================
if (WIN32)
  if ( NOT "${CMAKE_GENERATOR}" MATCHES "MinGW Makefiles" )
    add_custom_target (fixbundle ALL
      COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/fixbundle.win.cmake" )
  endif()
endif()

if (APPLE)
  add_custom_target (fixbundle ALL
    COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/fixbundle.cmake"
  )
endif()

# ========================================================================
# Bundle generation (Targets exist, now configure them)
# ========================================================================
if (WIN32)
  # prepare bundle generation cmake file and add a build target for it
  configure_file ("${OPENMESH_DIR}/cmake/fixbundle.cmake.win.in"
     "${CMAKE_BINARY_DIR}/fixbundle.win.cmake" @ONLY IMMEDIATE)

  if ( NOT "${CMAKE_GENERATOR}" MATCHES "MinGW Makefiles" )
    # let bundle generation depend on all targets
    add_dependencies (fixbundle QtViewer)
  endif()

endif()


if (APPLE)
  # prepare bundle generation cmake file and add a build target for it
  configure_file ("${OPENMESH_DIR}/cmake/fixbundle.cmake.in"
     "${CMAKE_BINARY_DIR}/fixbundle.cmake" @ONLY IMMEDIATE)

  # let bundle generation depend on all targets
  add_dependencies (fixbundle QtViewer) 

  # Required for Snow leopard, and the latest qt. Then the resources have to be copied
  if ( EXISTS "/opt/local/libexec/qt4-mac/lib/QtGui.framework/Versions/4/Resources/qt_menu.nib" )
      add_custom_command(TARGET OpenMesh POST_BUILD
         COMMAND ${CMAKE_COMMAND} -E copy_directory "/opt/local/libexec/qt4-mac/lib/QtGui.framework/Versions/4/Resources/qt_menu.nib" "${CMAKE_BINARY_DIR}/Build/Libraries/qt_menu.nib" 
        )
  endif ()

endif ()

# find needed packages for gui applications
find_package (OpenGL)
find_package (GLUT)

# For the apps, we need qt and opengl to build them
if (NOT QT4_FOUND)
  find_package (Qt4 COMPONENTS QtCore QtGui )

  set (QT_USE_QTOPENGL 1)

  include (${QT_USE_FILE})
endif ()

if ("${CMAKE_GENERATOR}" MATCHES "MinGW Makefiles")
  message(WARNING "GUI Apps are not build with mingw. (TODO)")
endif()

# check for OpenGL and GLUT as our required dependencies
if (QT4_FOUND AND OPENGL_FOUND AND GLUT_FOUND AND NOT "${CMAKE_GENERATOR}" MATCHES "MinGW Makefiles" )
  # Add ui apps as dependency before fixbundle 
  if ( WIN32 AND NOT "${CMAKE_GENERATOR}" MATCHES "MinGW Makefiles")
    # let bundle generation depend on all targets
    add_dependencies (fixbundle QtViewer)
  endif()

  # Add ui apps as dependency before fixbundle 
  if (APPLE)
    # let bundle generation depend on all targets
    add_dependencies (fixbundle QtViewer)
  endif()

  if (WIN32)
    FILE(GLOB files_install_app_dlls "${CMAKE_BINARY_DIR}/Build/*.dll" )
    INSTALL(FILES ${files_install_app_dlls} DESTINATION . )
  endif()

else () # Qt, Opengl or glut not found
  
  if (NOT QT4_FOUND)
    message ("QT 4 not found! Skipping some apps.")
  endif ()

  if (NOT OPENGL_FOUND)
    message ("OpengGL not found! Skipping some apps.")
  endif ()

  if (NOT GLUT_FOUND)
    message ("GLUT not found! Skipping some apps.")
  endif ()
  
endif ()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


set (targetName QtViewer)
set (OPENMESH_INCLUDE_DIR "${OPENMESH_DIR}/src")
set (OPENMESH_LIBRARY_DIR "${OPENMESH_DIR}/build/Build/lib")

include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/AddIncludes.cmake)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../src
  ${OPENMESH_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${QT_INCLUDE_DIR}
)

# source code directories
set (directories
  ${directories}
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../src
)

# collect all header and source files
acg_append_files (headers "*.h" ${directories})
acg_append_files (headers "*.hpp" ${directories})
acg_append_files (sources "*.c" ${directories})
acg_append_files (sources "*.cpp" ${directories})
acg_append_files (sources "*.cxx" ${directories})
acg_append_files (ui "*.ui" ${directories})

# remove template cc files from source file list
acg_drop_templates (sources)

# generate uic and moc targets
acg_qt4_autouic (uic_targets ${ui})
acg_qt4_automoc (moc_targets ${headers})


if (WIN32)
  acg_add_executable (${targetName} WIN32 ${uic_targets} ${sources} ${headers} ${moc_targets})
  # link to qtmain library to get WinMain function for a non terminal app
  target_link_libraries (${targetName} ${QT_QTMAIN_LIBRARY})
else ()
  acg_add_executable (${targetName} ${uic_targets} ${sources} ${headers} ${moc_targets})
endif ()

target_link_libraries (${targetName}
  ${QT_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${GLUT_LIBRARIES}
)

if (UNIX)
target_link_libraries (${targetName}
  debug ${OPENMESH_LIBRARY_DIR}/libOpenMeshCore${CMAKE_DEBUG_POSTFIX}.so
  debug ${OPENMESH_LIBRARY_DIR}/libOpenMeshTools${CMAKE_DEBUG_POSTFIX}.so
)
  
target_link_libraries (${targetName}
  ${OPENMESH_LIBRARY_DIR}/libOpenMeshCore.so
  ${OPENMESH_LIBRARY_DIR}/libOpenMeshTools.so
)
endif()

if (WIN32)
target_link_libraries (${targetName}
  debug ${OPENMESH_LIBRARY_DIR}/OpenMeshCore${CMAKE_DEBUG_POSTFIX}.lib
  debug ${OPENMESH_LIBRARY_DIR}/OpenMeshTools${CMAKE_DEBUG_POSTFIX}.lib
)
  
target_link_libraries (${targetName}
  optimized ${OPENMESH_LIBRARY_DIR}/OpenMeshCore.lib
  optimized ${OPENMESH_LIBRARY_DIR}/OpenMeshTools.lib
)
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/AddLibraries.cmake)
