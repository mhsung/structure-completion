##################################
#### Additional include files ####
##################################

if (UNIX)
set (LIBRARY_ROOT_PATH "~/lib" CACHE PATH "The directory where the all library files can be found.")

include_directories (
  ${LIBRARY_ROOT_PATH}/glog-0.3.3/build/include
  ${LIBRARY_ROOT_PATH}/Ipopt-3.12.1/build/include/coin
)
endif()

if (WIN32)
set (LIBRARY_ROOT_PATH "C:/project/lib" CACHE PATH "The directory where the all library files can be found.")

include_directories (
  #${LIBRARY_ROOT_PATH}/glog-0.3.3/src/windows
  ${LIBRARY_ROOT_PATH}/Ipopt-3.11.0/include/coin
)

# source code directories
if (MSVC)
#  file (GLOB SHADER_SOURCES
#   "${CMAKE_CURRENT_LIST_DIR}/../shader/*.h"
#	"${CMAKE_CURRENT_LIST_DIR}/../shader/*.cpp")
#  source_group("Shader Files" FILES ${SHADER_SOURCES})

  file (GLOB IPOPT_SOURCES
    "${CMAKE_CURRENT_LIST_DIR}/../ipopt/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../ipopt/*.cpp")
  source_group("IPOPT" FILES ${IPOPT_SOURCES})
  
  file (GLOB TRWS_SOURCES
    "${LIBRARY_ROOT_PATH}/TRW_S-v1.3/*.h"
	"${LIBRARY_ROOT_PATH}/TRW_S-v1.3/*.cpp")
  source_group("TRW_S" FILES ${TRWS_SOURCES})
  
  file (GLOB EIGENQP_SOURCES
    "${LIBRARY_ROOT_PATH}/wingsit_QP/*.h"
	"${LIBRARY_ROOT_PATH}/wingsit_QP/*.cpp")
  source_group("Eigen_QP" FILES ${EIGENQP_SOURCES})
endif()
endif()

include_directories (
  ${CMAKE_CURRENT_LIST_DIR}/../src
  ${CMAKE_CURRENT_LIST_DIR}/../ipopt
  #${CMAKE_CURRENT_LIST_DIR}/../shader
  ${LIBRARY_ROOT_PATH}/ann-1.1.2/include
  ${LIBRARY_ROOT_PATH}/eigen-3.2.2
  ${LIBRARY_ROOT_PATH}/gflags-2.1.1/build/include
  ${LIBRARY_ROOT_PATH}/glew-1.11.0/include
  ${LIBRARY_ROOT_PATH}/TRW_S-v1.3
  ${LIBRARY_ROOT_PATH}/wingsit_QP
  #${LIBRARY_ROOT_PATH}/ceres-solver-1.10.0/config
  #${LIBRARY_ROOT_PATH}/ceres-solver-1.10.0/include
)

set (directories
  ${CMAKE_CURRENT_LIST_DIR}/../ipopt
  #${CMAKE_CURRENT_LIST_DIR}/../shader
  ${LIBRARY_ROOT_PATH}/TRW_S-v1.3
  ${LIBRARY_ROOT_PATH}/wingsit_QP
)