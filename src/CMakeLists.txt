include(GNUInstallDirs)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

message(STATUS "BOOST_ROOT=${BOOST_ROOT}")
find_package(Boost COMPONENTS serialization REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LINK_DIR})

include_directories(${CMAKE_SOURCE_DIR}/include)
add_library(spmtrx SHARED spmcrs.cc spmutil.cc)

install(TARGETS spmtrx DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/spmtrx DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.h")
