if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(libigl)

export(TARGETS Eigen3_Eigen igl_core FILE libigl-config.cmake)