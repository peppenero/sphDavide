cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
        message(STATUS " The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

find_package(OpenCAL 1.0 COMPONENTS  OpenCAL  REQUIRED)
find_package(GLUT REQUIRED)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})



if(NOT MSVC)
  find_library(M_LIB m)
else()
  # not needed on windows
  set(M_LIB "")
endif()

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${OPENCAL_INCLUDE_DIR}
${OPENGL_INCLUDE_DIRS}
        ${GLUT_INCLUDE_DIR}

)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)

file(GLOB SRC
    "include/*.h"
    "source/*.cpp"
    "source/*.c"
    "*.c"
    "*.cpp"
    "*.h"
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_executable(cal_heat ${SRC})
target_link_libraries(cal_heat ${VTK_LIBRARIES})
target_link_libraries(cal_heat ${OPENCAL_LIBRARIES} ${M_LIB}  ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})
