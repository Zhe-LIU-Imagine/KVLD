SET(ProjectName kvld)

INCLUDE_DIRECTORIES(
../OrsaHomography/
../OrsaHomography/demo/libImage
)

IF(UNIX)
#support of ARRAY<T,N>
SET(CMAKE_CXX_FLAGS "-std=c++0x")
ENDIF(UNIX)

ADD_LIBRARY(kvld kvld.cpp kvld.h algorithm.cpp algorithm.h image.hpp matrix.cpp matrix.h match.h)



