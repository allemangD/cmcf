# Conformalized Mean Curvature Flow in Slicer

3D Slicer module for visualizing and analyzing Conformalized Mean Curvature Flow
on models.

# Building

Requires a local build of Slicer. https://slicer.readthedocs.io/en/latest/developer_guide/build_instructions/index.html

```shell
$ git clone https://github.com/allemangd/cmcf
$ cd cmcf
$ cmake -B ./build -DSlicer_DIR=$SLICER_INNER_BUILD
$ cmake --build ./build --parallel
$ ./build/SlicerWithSlicerCMCF
```
