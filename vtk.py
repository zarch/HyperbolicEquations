# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:17:37 2012

@author: pietro
"""
from __future__ import print_function

VTKTXT = """# vtk DataFile Version 3.1
{desc}
ASCII
DATASET UNSTRUCTURED_GRID

POINTS {lenpnts} FLOAT
{pnts}

CELLS {lentri} {celldim}
3 {tri}

CELL_TYPES {lentri}
{celltyp}

POINT_DATA {lenpnts}
{scalar}
"""

SCALARDATA = """SCALARS {label} FLOAT
LOOKUP_TABLE default
{values}
"""


def exportVTK(fname, points, triangles, pointsdata, description = ''):
    pnts = '\n'.join(['{0} {1} {2}'.format(x,y,z) for x, y, z in points])
    tri = '\n3 '.join(['{0} {1} {2}'.format(x,y,z) for x, y, z in triangles])
    celltypes = [5,] * len(triangles)
    celltyp = ' '.join([str(i) for i in celltypes])
    data = []
    for key, vals in pointsdata.items():
        values = '\n'.join(['{0}'.format(v) for v in vals])
        data.append(SCALARDATA.format(label = key, values = values ))
    scalar = '\n'.join(data)
    lentri = len(triangles)
    lenpnts = len(points)
    celldim = 4*len(triangles)
    filename = file(fname,'w+')
    print(VTKTXT.format(desc=description, lenpnts=lenpnts, pnts=pnts,
                  lentri=lentri, tri=tri,
                  celldim=celldim, celltyp=celltyp, scalar=scalar),
          file = filename)

p = [[0,0,0],[1,1,1],[1,0,1],[0,1,1]]
tri = [[0,1,2],[0,1,3]]
pdata = {'temperature' : [0,]*len(p),
         'pressure'    : [1,]*len(p)}
print('start')
exportVTK('prova.vtk', p, tri, pdata)
print('done!')

