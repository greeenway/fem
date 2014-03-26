import time
import numpy as np
import matplotlib.colors as mcolors

START_KEYWORD = '[BEGIN] '
END_KEYWORD   = '[END]   '

class Timer(object):
    """The Timer class times the durations of several statements """
    def __init__(self, name=None):
        self.name = name


    def __enter__(self):
        phrase = START_KEYWORD
        if self.name:
            phrase += self.name

        print(phrase)
        self.tstart = time.time()


    def __exit__(self, type, value, traceback):
        phrase = END_KEYWORD
        if self.name:
            phrase += self.name
        print('{0:30} ({1:5.3f}ms)'.format(phrase, (time.time() - self.tstart)*1000))

def write_matrix_to_file(matrix, filename):
    "Writes a numpy-matrix to a text file"
    f_hendl = open(filename, 'w')

    for row in matrix.tolist():
        i = 0
        for elem in row:
            f_hendl.write(str(elem))
            if i != len(row)-1:
               f_hendl.write(',\t') 
            i += 1

        f_hendl.write('\n')

    f_hendl.close()

def read_matrix_from_file(filename):
    "Reads a numpy-matrix from a text File"
    f_hendl = open(filename, 'r')

    lines = f_hendl.readlines()
    output = []

    for line in lines:
        output.append( [float(x) for x in line[:-1].split(',\t')])
    
    f_hendl.close()
    
    return np.array(output)

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


if __name__ == '__main__':
    with Timer('processing'):
        blup = 2
    with Timer('drawing'):
        print("blup")
    with Timer():
            pass



    A = np.random.random((4,3))
    print(A)
    write_matrix_to_file(A, 'data/test.txt')
    B = read_matrix_from_file('data/test.txt')
    print(B-A)
    C = read_matrix_from_file('data/nodes_to_coordinates.txt')
    print(C)

