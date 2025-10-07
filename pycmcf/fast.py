import asyncio
from subprocess import run
from tempfile import NamedTemporaryFile

import igl
import numpy as np
import scipy
import tqdm
import vtk
from scipy.sparse.linalg import spsolve
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk

# V: dense (v, 3)
# F: dense (f, 3)
# L: sparse
# M: sparse
# cotmatrix(V, F, L)
# massmatrix(V, F, ::MASSMATRIX_TYPE_BARYCENTRIC, M)
# V = solve(M - rate * L, M * V)

# r = vtk.vtkSTLReader()
# r.file_name = 'Stanford_Bunny.stl'

r = vtk.vtkPolyDataReader()
r.file_name = '/home/allem/src/cmcf/stx_stx_MNBCP116056-v02-1-6mo-20170324_t1w_bs_gray_surface_rsl_left_327680_native.vtk'

r = vtk.vtkCleanPolyData(input_connection=r.output_port)

r.Update()

O = np.reshape(
    r.output.GetPolys().GetConnectivityArray(),
    (-1,)
)[:-1]

F = np.reshape(
    r.output.GetPolys().GetConnectivityArray(),
    (-1, 3),
)

V = vtk_to_numpy(r.output.GetPoints().GetData()).copy()
V -= (center := np.mean(V, axis=0))
V /= (scale := np.mean(np.linalg.norm(V, axis=1)))

o = vtk.vtkPolyData()
o.DeepCopy(r.output)
o.GetPointData().RemoveArray(0)


async def flow(initial_rate=0.0001, acceleration=1.4, stages=25):
    global V

    r = initial_rate
    a = acceleration

    L = igl.cotmatrix(V, F)
    for stage in tqdm.tqdm(range(stages)):
        M = igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_BARYCENTRIC)
        V = scipy.sparse.linalg.spsolve(M - r * L, M * V, 'MMD_AT_PLUS_A')
        V /= np.mean(np.linalg.norm(V, axis=1))
        r *= a


async def show(V, *, reverse: bool = False, edges: bool = False):
    o.GetPoints().SetData(numpy_to_vtk(V * scale + center))

    with NamedTemporaryFile(suffix='.vtk', delete_on_close=False, delete=False) as f:
        f.close()
        print(f.name)
        w = vtk.vtkPolyDataWriter(file_name=f.name)

        if reverse:
            r = vtk.vtkReverseSense(input_data=o)
            w.input_connection = r.output_port
        else:
            w.input_data = o

        w.Update()

        OPTS = [
            '--filename=0',
            '--up=+Z',
            '--grid=0',
            '--axis=0',
            '--ambient-occlusion=1',
            '--tone-mapping=0',
            '--backface-type=hidden',
            '--camera-orthographic=1',
        ]

        if edges:
            OPTS += ['--edges=1']

        if reverse:
            OPTS += [
                '--camera-direction=-1,0,0',
                '--camera-position=300,0,0',
            ]
        else:
            OPTS += [
                '--camera-direction=1,0,0',
                '--camera-position=-300,0,0',
            ]

        run([ 'f3d', f.name, *OPTS ])


async def main():
    await flow(initial_rate=1, acceleration=10, stages=3)
    await show(V, edges=True)


asyncio.run(main())
