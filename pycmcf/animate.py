import asyncio
from concurrent.futures.thread import ThreadPoolExecutor
from pathlib import Path
from subprocess import run

import igl
import numpy as np
import scipy
import tqdm
import vtk
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

OUT = Path('out')
OUT.mkdir(exist_ok=True)


def write(stage):
    o.GetPoints().SetData(numpy_to_vtk(V * scale + center))

    w = vtk.vtkPolyDataWriter(
        file_name=OUT.joinpath(f'stage{stage:03}.vtk'),
        input_data=o,
    )
    w.Update()

    r = vtk.vtkReverseSense(
        input_data=o,
    )

    w = vtk.vtkPolyDataWriter(
        file_name=OUT.joinpath(f'stage{stage:03}rev.vtk'),
        input_connection=r.output_port,
    )
    w.Update()

    print(f'Wrote {stage:03}')


async def flow(stages=25):
    global V

    write(0)

    R = 0.0001
    A = 1.40

    L = igl.cotmatrix(V, F)

    for stage in tqdm.tqdm(range(stages)):
        M = igl.massmatrix(V, F, igl.MASSMATRIX_TYPE_BARYCENTRIC)
        V = scipy.sparse.linalg.spsolve(M - R * L, M * V, 'MMD_AT_PLUS_A')
        V /= np.mean(np.linalg.norm(V, axis=1))
        R *= A

        write(stage + 1)


async def animate():
    jobs = []

    OPTS = [
        '--filename=0',
        '--up=+Z',
        '--grid=0',
        '--axis=0',
        '--ambient-occlusion=1',
        '--tone-mapping=0',
        '--backface-type=hidden',
    ]

    with ThreadPoolExecutor(max_workers=4) as ex:
        for stage in range(26):
            jobs.append(asyncio.wrap_future(ex.submit(
                run,
                [
                    'f3d',
                    OUT.joinpath(f'stage{stage:03}.vtk'), '--output',
                    OUT.joinpath(f'stage{stage:03}.png'),
                    '--camera-direction=1,0,0',
                    '--camera-position=-300,0,0',
                    '--camera-orthographic=1',
                    *OPTS,
                ])))

            jobs.append(asyncio.wrap_future(ex.submit(run, [
                'f3d',
                OUT.joinpath(f'stage{stage:03}rev.vtk'), '--output',
                OUT.joinpath(f'stage{stage:03}-rev.png'),
                '--camera-direction=-1,0,0',
                '--camera-position=300,0,0',
                '--camera-orthographic=1',
                *OPTS,
            ])))

            jobs.append(asyncio.wrap_future(ex.submit(run, [
                'f3d',
                OUT.joinpath(f'stage{stage:03}.vtk'), '--output',
                OUT.joinpath(f'stage{stage:03}-edges.png'),
                '--camera-direction=1,0,0',
                '--camera-position=-300,0,0',
                '--camera-orthographic=1',
                '--edges=1',
                *OPTS,
            ])))

            jobs.append(asyncio.wrap_future(ex.submit(run, [
                'f3d',
                OUT.joinpath(f'stage{stage:03}rev.vtk'), '--output',
                OUT.joinpath(f'stage{stage:03}-edges-rev.png'),
                '--camera-direction=-1,0,0',
                '--camera-position=300,0,0',
                '--camera-orthographic=1',
                '--edges=1',
                *OPTS,
            ])))

    await asyncio.gather(*jobs)

    jobs = []
    for suf in ['', '-rev', '-edges', '-edges-rev']:
        jobs.append(asyncio.create_subprocess_exec(
            'ffmpeg', '-framerate', '5', '-i',
            f'out/stage%03d{suf}.png',
            '-loop', '1',
            '-crf', '19', f'stages{suf}.webm', '-y',
        ))

    await asyncio.gather(*jobs)


async def main():
    # await flow()
    await animate()


asyncio.run(main())
