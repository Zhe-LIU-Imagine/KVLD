CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(sift)

SET(SRC library.cpp library.h
        domain.cpp domain.h
        splines.cpp splines.h
        filter.cpp filter.h
        demo_lib_sift.cpp demo_lib_sift.h
        )

INCLUDE_DIRECTORIES(..)

ADD_LIBRARY(sift ${SRC})
IF(UNIX)
    SET_TARGET_PROPERTIES(sift PROPERTIES
                          COMPILE_FLAGS "-Wall -Wextra -Werror")
ENDIF(UNIX)

ADD_EXECUTABLE(sift_demo demo_sift.cpp)
TARGET_LINK_LIBRARIES(sift_demo sift image Numerics)

INSTALL(FILES CMakeLists.txt ${SRC} demo_sift.cpp DESTINATION src/extras/sift)
