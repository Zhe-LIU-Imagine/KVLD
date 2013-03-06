SET(ProjectName KVLD_DEMO)
PROJECT(${ProjectName})

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(OpenMP) 
if (OPENMP_FOUND) 
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

#Detect c++11 support on Unix system
if (UNIX)
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
      execute_process(
          COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
      if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
          message(STATUS "C++11 activated.")
          add_definitions("-std=gnu++11")
        elseif (GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
          message(WARNING "C++0x activated. If you get any errors update to a compiler which fully supports C++11")
          add_definitions("-std=gnu++0x")
        endif ()
      endif ()
  elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  else ()
      message(FATAL_ERROR "Your C++ compiler does not support C++11.")
  endif ()
endif (UNIX)

# If you don't want to use openCV  and only need KVLD, 
# please commande the next line, 
# don't forget to also commande ${OpenCV_LIBS} further.
FIND_PACKAGE( OpenCV REQUIRED )


#find Orsa and KVLD librarys
ADD_DEFINITIONS(-DSOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
ADD_SUBDIRECTORY(OrsaHomography)
ADD_SUBDIRECTORY(kvld)

INCLUDE_DIRECTORIES(
./OrsaHomography/
./OrsaHomography/demo/libImage
$ENV{OPENCV_DIR}/install/include/opencv
)


ADD_EXECUTABLE(KVLD_Deformable deformable.cpp convert.h convert.cpp )
TARGET_LINK_LIBRARIES(KVLD_Deformable kvld image sift Orsa
${OpenCV_LIBS})

ADD_EXECUTABLE(KVLD_Calibration calibration.cpp convert.h convert.cpp)
TARGET_LINK_LIBRARIES(KVLD_Calibration kvld image sift Orsa
${OpenCV_LIBS})

ADD_EXECUTABLE(KVLD_Interface interface.cpp convert.h convert.cpp)
TARGET_LINK_LIBRARIES(KVLD_Interface kvld image sift Orsa
${OpenCV_LIBS})


