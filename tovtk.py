import sys, os
from glob import glob
import numpy as np
from myio import Reader
import pickle

def read(filename):
    reader = Reader()
    Q1, Q2, Q3 = reader.read(filename)
    q1 = Q1.reshape([nx,ny], order='F')
    q2 = Q2.reshape([nx,ny], order='F')
    q3 = Q3.reshape([nx,ny], order='F')
    I = q1**2+q2**2+q3**2
    s1 = q2*q3
    s2 = -q1*q3
    s3 = 0.0*q1
    S = s1**2 + s2**2
    return q1, q2, q3, s1, s2, s3, I, S

def tobin(state,varname,base_name,path='./'):
    filename = varname+'_'+base_name
    output_file  = open(os.path.join(path,filename),'wb')
    state.tofile(output_file)
    output_file.close()

def tovtk(q1,q2,q3,s1,s2,s3,I,S,base_name,path):
    scalars = [("Q1", q1),
               ("Q2", q2),
               ("Q3", q3),
               ("S1", s1),
               ("S2", s2),
               ("S3", s3),
               ("I",   I),
               ("S",   S)]
    vectors = []
    title = 'VTK Data'
    filename = os.path.join(path,base_name+'.vtk')
    print filename
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

if __name__ == "__main__":
    import sys
    path = sys.argv[1]
    # nx = int(sys.argv[2])
    # ny = int(sys.argv[3])
    # print nx,ny    
    # # nx, ny = 1000, 221


    # get the pickle file and cell number

    pkl_file_name = glob(os.path.join(path,'*.pkl'))
    pkl_file = open(pkl_file_name[0],'rb')
    pkl = pickle.load(pkl_file)

    nx = pkl.get('nx')
    ny = pkl.get('ny')
    print nx,ny
    
    coordinates = [np.linspace(0,1,nx),
                   np.linspace(0,1,ny),
                   np.ones(1),
                   ]
    dimensions = (nx, ny, 1) 

    file_name = glob(os.path.join(path,'step*.dat'))
    file_name.sort()
    for n,filename in enumerate(file_name):
        m = filename.split('step')[1].split('.')[0].zfill(7)
        print m
        q1,q2,q3,s1,s2,s3,I,S = read(filename)
        # vtkname = filename.replace('.dat', '.vtk')
        # base_name = str(m).zfill(4)+'.databin'
        # tobin(s1,'s1',base_name,path)
        # tobin(s2,'s2',base_name,path)
        # tobin(q1,'q1',base_name,path)
        # tobin(q2,'q2',base_name,path)
        # tobin(q3,'q3',base_name,path)
        # tobin(S,'S',base_name,path)
        # tobin(I,'I',base_name,path)

        base_name = str(m).zfill(4)
        tovtk(q1,q2,q3,s1,s2,s3,I,S,base_name,path)
