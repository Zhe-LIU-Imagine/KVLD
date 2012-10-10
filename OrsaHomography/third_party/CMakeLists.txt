#---------------------
#Third party libraries
#---------------------

ADD_SUBDIRECTORY(CppUnitLite)
ADD_SUBDIRECTORY(svd)

SET(INCLUDE_INSTALL_DIR src/third_party)
#ADD_SUBDIRECTORY(Eigen)

ADD_SUBDIRECTORY(zlib)
ADD_SUBDIRECTORY(png)
ADD_SUBDIRECTORY(jpeg)

INSTALL(FILES CMakeLists.txt DESTINATION src/third_party)
