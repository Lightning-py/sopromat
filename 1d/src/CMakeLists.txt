cmake_minimum_required(VERSION 3.10)
project(MKE2d VERSION 0.1.0 LANGUAGES CXX)

# set(SOURCE_DIR src)

set(HEADERS headers)
set(SOURCES sources)

set(SOURCE_FILES ${SOURCES}/dot.cpp ${SOURCES}/element.cpp ${SOURCES}/force.cpp)

include_directories(${HEADERS})


# set(GLUTILS ${HEADERS}/glUtils.hpp ${SOURCES}/glUtils.cpp)
# set(RANDOMUTILS ${HEADERS}/randomUtils.hpp ${SOURCES}/randomUtils.cpp)

# add_library(LIB_UTILS STATIC ${GLUTILS} ${RANDOMUTILS})

# add_link_options(-lglut -lGL -lGLU -lGLEW)
add_executable(MKE2d main.cpp ${SOURCE_FILES})

# target_link_libraries(FlappyBird LIB_UTILS)