cmake_minimum_required(VERSION 3.8)
project(spin_dynamics)

#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS -std=c++14")

file( GLOB LIB_SOURCES lib/*.cpp )
file( GLOB LIB_HEADERS lib/*.h )
add_library( Lib   ${LIB_SOURCES} ${LIB_HEADERS})
#set_target_properties(Lib PROPERTIES LINKER_LANGUAGE CXX)

#Eigen
set(EIGEN_INCLUDE_DIR "" CACHE PATH "Include directory for Eigen 3.")
include_directories(SYSTEM ${EIGEN_INCLUDE_DIR} )

#Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(Boost_USE_STATIC_LIBS OFF)
set (Boost_INCLUDE_DIR "" CACHE PATH "Include directory for Boost.")
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

configure_file(grid/make_grid.m grid/make_grid.m COPYONLY)
file(COPY "grid/" DESTINATION "grid/" FILES_MATCHING PATTERN "*.txt")


add_executable(Echo_Doublet Echo_Doublet.cpp  )
target_link_libraries( Echo_Doublet Lib )

add_executable(EDFS_Doublet EDFS_Doublet.cpp  )
target_link_libraries( EDFS_Doublet  Lib )

add_executable(Echo_Triplet Echo_Triplet.cpp  )
target_link_libraries( Echo_Triplet Lib )

add_executable(EDFS_Nitroxide EDFS_Nitroxide.cpp  )
target_link_libraries( EDFS_Nitroxide Lib )

add_executable(EDFS_Doublet_Triplet EDFS_Doublet_Triplet.cpp)
target_link_libraries( EDFS_Doublet_Triplet Lib )

add_executable(Echo_Nitroxide Echo_Nitroxide.cpp  )
target_link_libraries( Echo_Nitroxide Lib )

add_executable(EDFS_Triplet EDFS_Triplet.cpp  )
target_link_libraries( EDFS_Triplet Lib )

add_executable(LaserIMD LaserIMD.cpp  )
target_link_libraries( LaserIMD Lib )


