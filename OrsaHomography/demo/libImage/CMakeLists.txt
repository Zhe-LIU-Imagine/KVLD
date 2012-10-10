
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

SET(SRC image.hpp
		sample.hpp
		pixelTypes.hpp
		image_converter.hpp
		image_io.hpp
		image_io.cpp
		image_drawing.hpp
		image_concat.hpp)

ADD_LIBRARY(image ${SRC})

TARGET_LINK_LIBRARIES(image png jpeg)

IF(UNIX)
    SET_TARGET_PROPERTIES(image PROPERTIES
                          COMPILE_FLAGS "-Wall -Wextra -Werror")
ENDIF(UNIX)

UNIT_TEST(image "image")
UNIT_TEST(image_drawing "image")
UNIT_TEST(image_io "image")

INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/demo/libImage)
INSTALL(FILES image_test.cpp image_drawing_test.cpp image_io_test.cpp
        DESTINATION src/demo/libImage)
SET(IMG test_image/two_pixels_color.jpg test_image/two_pixels_monochrome.jpg test_image/two_pixels.ppm  test_image/two_pixels_color.png test_image/two_pixels_monochrome.png test_image/two_pixels_gray.pgm test_image/two_pixels.pgm)

INSTALL(FILES ${IMG} DESTINATION src/demo/libImage/image_test)
