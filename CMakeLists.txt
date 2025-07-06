cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(SlicerCMCF)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/allemangD/cmcf")
set(EXTENSION_CONTRIBUTORS "David Allemang (University of North Carolina at Chapel Hill)")
set(EXTENSION_DESCRIPTION "Conformalized mean curvature flow (CMCF) implementation and visualizations in Slicer.")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/SlicerCMCF.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/SlicerCMCF/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CMCFlib)
add_subdirectory(CMCF)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
