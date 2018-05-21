project (griddata)


cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
find_package(CGAL REQUIRED)
find_package(OpenMP)
include( ${CGAL_USE_FILE} )
set(CMAKE_INCLUDE_CURRENT_DIR ON)
#find_package(HDF5 REQUIRED)

set (LIB_TYPE STATIC) # or SHARED
string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)
find_package (HDF5 NAMES hdf5 COMPONENTS C CXX ${SEARCH_TYPE})
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_${LIB_TYPE}_LIBRARY} ${HDF5_CXX_${LIB_TYPE}_LIBRARY})

file(GLOB HDRS	"*.h")
file(GLOB SRCS	"*.cpp")


add_executable(griddata ${SRCS} ${HDRS})
#TARGET_C_PROPERTIES (griddata ${LIB_TYPE} " " " ")
message(${LINK_LIBS})
target_link_libraries (griddata ${LINK_LIBS})

