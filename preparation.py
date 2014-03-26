import numpy as np
import math
from math import floor as floor




def create_local_to_global(info):
    """
    Creates the local_to_global table. Maps local to global coordinates for
    each element.
    """
    n = info.elements_per_line
    N = info.number_of_elements

    local_to_global = np.zeros(shape=(N, info.nodes_per_element)) #one colum less

    for k in range(1,N+1):
        i = floor((k-1.)/n) + 1

        if k-1 == 0:
            i = 1

        local_to_global[k-1,0] = (k-1+i);
        local_to_global[k-1,1] = (k-0+i);
        local_to_global[k-1,2] = (k+1+i+n);
        local_to_global[k-1,3] = (k+i+n);

    return local_to_global


def create_nodes_to_coordinates(info):
    """
    Create nodes_to_coordinates table. Maps the global nodes to their 
    (x,y)-coordinates.
    """
    N_nodes = info.number_of_nodes
    n = info.elements_per_line

    nodes_to_coordinates = np.zeros(shape=(N_nodes, 2))

    for k in range(0,N_nodes):
        nodes_to_coordinates[k,0] = (k % (n+1)) * 2; #x
        nodes_to_coordinates[k,1] = floor((k/(n+1))) * 2; #y

    return nodes_to_coordinates


def create_boundary_table(info):
    """
    Create boundary_table. Specifies the boundary conditions for
    each global node.
    """

    N_nodes = info.number_of_nodes
    n = info.elements_per_line
    n1 = n+1

    #boundary_table = np.zeros(shape=(2, N_nodes))

    #matrix form
    set1 = np.zeros(shape=(n1,n1))
    value1 = np.zeros(shape=(n1,n1))

    for i in range(0,n1):
        for j in range(0,n1):
            # (i,j) "coordinates"
            if  i == 0 or i == n or j == 0 or j == n:
                #set border to zero
                set1[i,j] = 1
                value1[i,j] = 0



            # if i == 5 and j > 4 and j <= 8:
            #     set1[i,j] = 1
            #     value1[i,j] = 2

            # if i == 7 and j > 4 and j <= 8:
            #     set1[i,j] = 1
            #     value1[i,j] = -2


            

    #convert matrix form to "lines"
    l1 = set1.reshape(1,N_nodes)   #first line dirichlet set
    l2 = value1.reshape(1,N_nodes) #second line dirichlet value

    return np.vstack((l1,l2))




#execute run function when executed
if __name__ == '__main__':
    import fem
    fem.run()
































