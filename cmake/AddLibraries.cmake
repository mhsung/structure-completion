##################################
#### Additional library files ####
##################################

if (UNIX)
target_link_libraries (${targetName}
  ${LIBRARY_ROOT_PATH}/ann-1.1.2/lib/libANN.a
  ${LIBRARY_ROOT_PATH}/glew-1.11.0/lib/libGLEW.a
  ${LIBRARY_ROOT_PATH}/glew-1.11.0/lib/libGLEW.so
)

target_link_libraries (${targetName}
  ${LIBRARY_ROOT_PATH}/gflags-2.1.1/build/lib/libgflags.a
  ${LIBRARY_ROOT_PATH}/glog-0.3.3/build/lib/libglog.a
  ${LIBRARY_ROOT_PATH}/glog-0.3.3/build/lib/libglog.so
  ${LIBRARY_ROOT_PATH}/ceres-solver-1.10.0/build/lib/libceres.a
)
endif()

if (WIN32)
target_link_libraries (${targetName}
  ${LIBRARY_ROOT_PATH}/ann-1.1.2/build/bin/ANN.lib
  ${LIBRARY_ROOT_PATH}/glew-1.11.0/lib/Release/x64/glew32.lib
)

target_link_libraries (${targetName}
  debug ${LIBRARY_ROOT_PATH}/gflags-2.1.1/build/lib/Debug/gflags.lib
  debug ${LIBRARY_ROOT_PATH}/glog-0.3.3/x64/Debug/libglog.lib
  debug ${LIBRARY_ROOT_PATH}/ceres-solver-1.10.0/build/lib/Debug/ceres-debug.lib
)
  
target_link_libraries (${targetName}
  optimized ${LIBRARY_ROOT_PATH}/gflags-2.1.1/build/lib/Release/gflags.lib
  optimized ${LIBRARY_ROOT_PATH}/glog-0.3.3/x64/Release/libglog.lib
  optimized ${LIBRARY_ROOT_PATH}/ceres-solver-1.10.0/build/lib/Release/ceres.lib
)
endif()
