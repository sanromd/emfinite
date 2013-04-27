import sys, os
from glob import glob
import numpy as np
from io import Reader

def read(filename):
    reader = Reader()
    Q1, Q2, Q3 = reader.read(filename)
    return Q1, Q2, Q3
    q1 = Q1.reshape([nx,ny], order='F')
    q2 = Q2.reshape([nx,ny], order='F')
    q3 = Q3.reshape([nx,ny], order='F')

nx, ny = 35, 35

coordinates = [np.linspace(0,1,nx),
               np.linspace(0,1,ny),
               np.ones(1),
               ]
dimensions = (nx, ny, 1) 

def tovtk(q1,q2,q3,filename):
    scalars = [("Q1", q1),
               ("Q2", q2),
               ("Q3", q3)]
    vectors = []
    title = 'VTK Data'

    fh = open(filename, 'wb')
    fh_write = lambda s: fh.write(s.encode('ascii'))

    header = '# vtk DataFile Version %d.%d'
    version = (2, 0)
    fh_write(header % version)
    fh_write('\n')
    title = title
    fh_write(title[:255])
    fh_write('\n')

    format = 'BINARY'
    fh_write(format)
    fh_write('\n')

    dataset_type = 'RECTILINEAR_GRID'
    fh_write('DATASET %s' % dataset_type);
    fh_write('\n')
    fh_write('DIMENSIONS %d %d %d' % dimensions)
    fh_write('\n')
    for X, array in zip("XYZ", coordinates):
        label = X+'_COORDINATES'
        fh_write('%s %s %s' % (label, len(array), 'double'))
        fh_write('\n')
        array.astype('>d').tofile(fh)
        fh_write('\n')

    data_type = 'POINT_DATA'
    fh_write('%s %d' % (data_type, np.prod(dimensions)))
    fh_write('\n')

    for i, (name, array) in enumerate(scalars):
        attr_type = 'SCALARS'
        attr_name = name or (attr_type.lower() + str(i))
        attr_name = attr_name.replace(' ', '_')
        fh_write('%s %s %s' %(attr_type, attr_name, 'double'))
        fh_write('\n')
        lookup_table = 'default'
        lookup_table = lookup_table.replace(' ', '_')
        fh_write('LOOKUP_TABLE %s' % lookup_table)
        fh_write('\n')
        array.astype('>d').tofile(fh)
        fh_write('\n')

    for i, (name, array) in enumerate(vectors):
        attr_type = 'VECTORS'
        attr_name = name or (attr_type.lower() + str(i))
        attr_name = attr_name.replace(' ', '_')
        fh_write('%s %s %s' %(attr_type, attr_name, 'double'))
        fh_write('\n')
        array.astype('>d').tofile(fh)
        fh_write('\n')

    fh.flush()
    fh.close()


for filename in glob("step*.dat"):
    q1,q2,q3 = read(filename)
    vtkname = filename.replace('.dat', '.vtk')
    tovtk(q1,q2,q3,vtkname)
