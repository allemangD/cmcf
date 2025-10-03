import slicer.util
import vtk
import numpy as np

for side in ['left', 'right']:
    orig = slicer.util.getNode(f'*{side}*native')

    orig_filt = vtk.vtkCellSizeFilter()
    orig_filt.SetInputData(orig.GetPolyData())
    orig_filt.Update()
    orig_area = vtk.vtkDoubleArray()
    orig_area.DeepCopy(orig_filt.GetOutput().GetCellData().GetArray('Area'))

    post = slicer.util.getNode(f'*{side}*CMCF')

    post_filt = vtk.vtkCellSizeFilter()
    post_filt.SetInputData(post.GetPolyData())
    post_filt.Update()
    post_area = vtk.vtkDoubleArray()
    post_area.DeepCopy(post_filt.GetOutput().GetCellData().GetArray('Area'))

    post_area_array = np.asarray(post_area)
    post_area_array /= np.asarray(orig_area)

    post.GetPolyData().GetCellData().AddArray(post_area)

    copy = slicer.mrmlScene.CopyNode(post)
    copy.SetName(f'final cmcf ratios - {side}')
