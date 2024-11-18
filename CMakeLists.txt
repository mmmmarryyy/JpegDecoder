file(GLOB_RECURSE SOLUTION_SRC CONFIGURE_DEPENDS "decoder/*.cpp")
add_shad_library(decoder ${SOLUTION_SRC})

add_catch(test_baseline tests/test_commons.cpp tests/test_baseline.cpp)
add_catch(test_progressive tests/test_commons.cpp tests/test_progressive.cpp)

set(LIBRARIES decoder ${FFTW_LIBRARIES} ${PNG_LIBRARY} ${JPEG_LIBRARIES})
target_link_libraries(test_baseline PRIVATE ${LIBRARIES})
target_link_libraries(test_progressive PRIVATE ${LIBRARIES})
