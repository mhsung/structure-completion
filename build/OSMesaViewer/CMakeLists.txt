cmake_minimum_required (VERSION 2.8)
set (targetName OSMesaViewer)

project (${targetName})

include (${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/FindOSMesa.cmake)
if (NOT OSMESA_FOUND)
  message (FATAL_ERROR "OSMesa not found!")
endif ()

# import common lists
include (${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/CMakeListsCommon.cmake)

target_include_directories (${targetName} PUBLIC
  ${OSMESA_INCLUDE_DIR}
)

target_link_libraries (${targetName}
  ${OSMESA_LIBRARY}
)

