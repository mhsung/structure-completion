cmake_minimum_required (VERSION 2.8)
set (target_name OSMesaViewer)

project (${target_name})

include (${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/FindOSMesa.cmake)
if (NOT OSMESA_FOUND)
  message (FATAL_ERROR "OSMesa not found!")
endif ()

# import common lists
include (${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/CMakeListsCommon.cmake)

target_include_directories (${target_name} PUBLIC
  ${OSMESA_INCLUDE_DIR}
)

target_link_libraries (${target_name}
  ${OSMESA_LIBRARY}
)

