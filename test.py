#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib
import itertools
import math

import numeric


phi_lin_tri = [
    lambda xi1, xi2: 1-xi1-xi2,
    lambda xi1, xi2: xi1,
    lambda xi1, xi2: xi2,
]

def _read_edges(filename):
    """
    Converts 'e'-csv dump to ?

    """
    f = open(filename, 'r')

    e1 = np.array([int(x) for x in f.readline()[:-1].split(',')]) #1. global node
    e2 = np.array([int(x) for x in f.readline()[:-1].split(',')]) #2. global node
    f.readline()
    f.readline()
    e5 = np.array([int(x) for x in f.readline()[:-1].split(',')]) # edge number

    f.close()

    edge2global = np.vstack((e1,e2)).transpose()-1 #0 index
    edge2number = e5.transpose()

    return edge2global, edge2number

def create_edges(filename, local_to_global):
    edge2global, edge2number = _read_edges(filename)

    edge_to_local = np.zeros(shape=edge2global.shape)
    tmp_row = np.zeros(shape=(1,edge_to_local.shape[1]))

    for i in range(0, edge_to_local.shape[0]):
        for j in range(0, local_to_global.shape[0]):
            found = 0
            for k in range(0, 3):
                for m in range(edge_to_local.shape[1]):
                    if local_to_global[j,k] == edge2global[i,m]:
                        tmp_row[0,m] = k
                        found += 1
                        #print('found!')
            if found == 2:
                edge_to_local[i,:] = tmp_row
                #print('found row!')
                break

    
    e2l = edge_to_local
    e2g = edge2global
    #l2g = local_to_global
    #print(l2g)
    #print(edge_to_local.shape[0])
    
    ##debuging
    # for i in range(0,edge_to_local.shape[0]):
    #     print("({0:d},{1:d}) -> ({2:d},{3:d})".format(int(e2l[i,0]), 
    #         int(e2l[i,1]), int(e2g[i,0]), int(e2g[i,1])))



    
    return edge2global, edge2number, edge_to_local


def _read_elements(filename):
    """
    Converts 't'-csv dump to a Matrix T

    """
    f = open(filename, 'r')

    t1 = np.array([int(x) for x in f.readline()[:-1].split(',')])
    t2 = np.array([int(x) for x in f.readline()[:-1].split(',')])
    t3 = np.array([int(x) for x in f.readline()[:-1].split(',')])

    f.close()

    T = np.vstack((t1,t2,t3))-1 #0 index

    return T

def create_local_to_global_tri(filename):
    return _read_elements(filename).transpose()


#

def create_nodes_to_coordinates_tri(filename):
    return _read_coordinates(filename).transpose()



def _read_coordinates(filename):
    """
    Converts 'p'-csv dump to a Matrix P

    """
    f = open(filename, 'r')

    p1 = np.array([float(x) for x in f.readline()[:-1].split(',')])
    p2 = np.array([float(x) for x in f.readline()[:-1].split(',')])

    f.close()

    P = np.vstack((p1,p2))

    return P

#borders...


def _ref_to_global(x1, x2, x3, y1, y2, y3, xi):
    xi1 = xi[0]
    xi2 = xi[1]

    #x,y = J_r * xi  + x1
    x = (x2 - x1) * xi1 + (x3 - x1) * xi2 + x1;
    y = (y2 - y1) * xi1 + (y3 - y1) * xi2 + y1;

    return (x,y)



def _calculate_K_r(x1, x2, x3, y1, y2, y3, eps):
    """
    _calculate_K_r
    """
    K_r = np.zeros(shape=(3,3))
    #gradients do not depend on x,y in this case
    #therefore K_r = J^-T*grad_i*J^-T*grad_j * abs(det_J) * int ( eps(...) )

    detJ = x2 * y3 - x2 * y1 - x1 * y3 - x3 * y2 + x3 * y1 + x1 * y2;

    K_r[0,0] = (y3 * y3 - 2 * y2 * y3 + y2 * y2 + x2 * x2 - 2 * x2 * x3 + x3 * x3)
    K_r[0,1] = -(y3 * y3 - y3 * y1 - y2 * y3 + y2 * y1 - x2 * x3 + x2 * x1 + x3 * x3 - x3 * x1)
    K_r[0,2] = (y2 * y3 - y2 * y2 - y3 * y1 + y2 * y1 + x2 * x3 - x2 * x2 - x3 * x1 + x2 * x1)
    K_r[1,0] = -(y3 * y3 - y3 * y1 - y2 * y3 + y2 * y1 - x2 * x3 + x2 * x1 + x3 * x3 - x3 * x1)
    K_r[1,1] = (y3 * y3 - 2 * y3 * y1 + y1 * y1 + x3 * x3 - 2 * x3 * x1 + x1 * x1)
    K_r[1,2] = -(y2 * y3 - y2 * y1 - y3 * y1 + y1 * y1 + x2 * x3 - x2 * x1 - x3 * x1 + x1 * x1)
    K_r[2,0] = (y2 * y3 - y2 * y2 - y3 * y1 + y2 * y1 + x2 * x3 - x2 * x2 - x3 * x1 + x2 * x1)
    K_r[2,1] = -(y2 * y3 - y2 * y1 - y3 * y1 + y1 * y1 + x2 * x3 - x2 * x1 - x3 * x1 + x1 * x1)
    K_r[2,2] = (y2 * y2 - 2 * y2 * y1 + y1 * y1 + x2 * x2 - 2 * x2 * x1 + x1 * x1)

    epsIntegral = numeric.gauss_triangle_6(lambda x,y: eps(_ref_to_global(x1, x2, x3, y1, y2, y3, (x,y))))
    
    K_r = K_r / detJ**2 * abs(detJ) * epsIntegral

    return K_r


def _calculate_f_r(x1, x2, x3, y1, y2, y3, rho):
    """
    _calculate_f_r elementwise "lastvektoren"
    """
    f_r = np.zeros(shape=(3,1))
    phi = phi_lin_tri
    adetJ = abs(x2 * y3 - x2 * y1 - x1 * y3 - x3 * y2 + x3 * y1 + x1 * y2)

    for i in range(0,3):
        integrand = lambda xi1, xi2: rho(_ref_to_global(x1, x2, x3, y1, y2, y3, (xi1,xi2))) * phi[i](xi1,xi2) * adetJ

        f_r[i,0] = numeric.gauss_triangle_6(integrand)

    return f_r



def assemble_equations(local_to_global, nodes_to_coordinates, borders, eps, rho, neumann):
    """
    Assembles the final system of linear equations to be solved.
    Outputs a Matrix K_h and the right hand side vector f_h
    """

    N_nodes = nodes_to_coordinates.shape[0]
    N_fe = local_to_global.shape[0]

    K_h = np.zeros(shape=(N_nodes, N_nodes))
    f_h = np.zeros(shape=(N_nodes,1))

    for r in range(0, N_fe):
        
        #K_r, f_r for each element
        i1 = local_to_global[r, 0]
        i2 = local_to_global[r, 1]
        i3 = local_to_global[r, 2]
        x1 = nodes_to_coordinates[i1, 0] 
        x2 = nodes_to_coordinates[i2, 0] 
        x3 = nodes_to_coordinates[i3, 0] 
        y1 = nodes_to_coordinates[i1, 1] 
        y2 = nodes_to_coordinates[i2, 1] 
        y3 = nodes_to_coordinates[i3, 1] 

        K_r = _calculate_K_r(x1, x2, x3, y1, y2, y3, eps)
        f_r = _calculate_f_r(x1, x2, x3, y1, y2, y3, rho) 

        for alpha in range(0,3):
            i = local_to_global[r,alpha]
            
            f_h[i,0] = f_h[i,0] + f_r[alpha,0]
            
            for beta in range(0,3):
                j = local_to_global[r,beta]

                K_h[i,j] = K_h[i,j] + K_r[alpha,beta] 

    #neumann boundarys
    #neumann = {"e2g": edge_to_global, "e2n": edge_to_number, "e2l": edge_to_local, "funcs":func_dict} 
    edge_to_global = neumann["e2g"]
    edge_to_number = neumann["e2n"]
    edge_to_local = neumann["e2l"]
    func_dict = neumann["funcs"]

    Ne = edge_to_local.shape[1]

    for e in range(0, edge_to_global.shape[0]):
        #only if e in neumann boundary TODO
        #print(func_dict.keys())
        
        if edge_to_number[e] in func_dict.keys():

            e_number = edge_to_number[e]

            #create fe2
            fe2 = np.zeros(shape=(Ne,1))
            
            xe_1 = nodes_to_coordinates[edge_to_global[e, 0], 0]
            xe_Ne = nodes_to_coordinates[edge_to_global[e, Ne-1], 0]
            ye_1 = nodes_to_coordinates[edge_to_global[e, 0], 1]
            ye_Ne = nodes_to_coordinates[edge_to_global[e, Ne-1], 1]

            root =  math.sqrt( (xe_Ne-xe_1)**2 + (ye_Ne-ye_1)**2 )

            #x = (xe_Ne - xe_1) + xi1 + xe_1
            #y = (ye_Ne - ye_1) + xi2 + ye_1
            x_e = lambda xi: (xe_Ne - xe_1) + xi + xe_1
            y_e = lambda xi: (ye_Ne - ye_1) + xi + ye_1

            #phi_lin_tri

            #int func(global_to_ref(xi)) phi_(p,e) root ds

            l1 = int(edge_to_local[e,0])
            l2 = int(edge_to_local[e,1])

            
            #do this with for loop
            fe2[0,0] = numeric.gauss_1D_6( lambda xi: func_dict[e_number](x_e(xi), y_e(xi)) * phi_lin_tri[l1](x_e(xi), y_e(xi)) * root )
            fe2[1,0] = numeric.gauss_1D_6( lambda xi: func_dict[e_number](x_e(xi), y_e(xi)) * phi_lin_tri[l2](x_e(xi), y_e(xi)) * root )

            for p in range(0, Ne):
                i = edge_to_global[e,p]
                f_h[i,0] = f_h[i,0] + fe2[p,0]



    #dirichlet boundarys
    #borders2 = [b[0] for b in borders]
    #borders3 = list(itertools.chain(*borders2))

    values = np.zeros(shape=(N_nodes,1))
    is_set = np.zeros(shape=(N_nodes,1))
    
    for nodes, value in borders:
        for n in nodes:
           is_set[n] = 1
           values[n] = value

    b_tmp = values * is_set

    for i in range(0, N_nodes):
        for j in range(0, N_nodes):
            if is_set[i] == 0:
                f_h[i] = f_h[i] - K_h[i,j] * b_tmp[j]

    for i in range(0,N_nodes):
        for j in range(0,N_nodes):
            if is_set[i] == 1 or is_set[j] == 1:
                if i == j:
                    K_h[i,j] = 1
                else:
                    K_h[i,j] = 0

        if is_set[i] == 1:
            f_h[i] = values[i]

    return K_h, f_h

def contour_plot(local_to_global, nodes_to_coordinates, u_h):


    x = nodes_to_coordinates[:,0]#.transpose()
    y = nodes_to_coordinates[:,1]#.transpose()
    triangulation = tri.Triangulation(x, y, local_to_global)

    N_fe = local_to_global.shape[0]

    u_max = u_h.max()
    u_min = u_h.min()

    umin2 = u_min-0.05*abs(u_max-u_min)
    umax2 = u_max+0.05*abs(u_max-u_min)

    levels = np.arange(u_min-0.05*abs(u_max-u_min), u_max+0.05*abs(u_max-u_min), abs(u_max-u_min)/25.)

    #levels = np.arange(-1.1, 1.1, 0.1)

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.8])
    ax.set_xlim((-1,1)) #fix constants
    ax.set_ylim((-1,1))

    for r in range(0,N_fe):
        i1 = local_to_global[r, 0]
        i2 = local_to_global[r, 1]
        i3 = local_to_global[r, 2]

        x1 = nodes_to_coordinates[i1, 0] 
        x2 = nodes_to_coordinates[i2, 0] 
        x3 = nodes_to_coordinates[i3, 0] 
        y1 = nodes_to_coordinates[i1, 1] 
        y2 = nodes_to_coordinates[i2, 1] 
        y3 = nodes_to_coordinates[i3, 1]

        t = tri.Triangulation(x[[i1, i2, i3]], y[[i1, i2, i3]])
        r = tri.UniformTriRefiner(t)

        # zX = t.x.copy()
        # zX[0] = u_h[i1]
        # zX[1] = u_h[i2]
        # zX[2] = u_h[i3]

        #zX = np.array([u_h[i1],u_h[i2],u_h[i3]])
        
        #print(t.x)
        #print(zX)
        refi_t = r.refine_triangulation(subdiv=2)

        rx = refi_t.x.copy()
        ry = refi_t.y.copy()


        detJ = x2 * y3 - x2 * y1 - x1 * y3 - x3 * y2 + x3 * y1 + x1 * y2

        rx2 = rx.copy()
        ry2 = ry.copy()

        for i in range(0,rx.shape[0]): #buggy?


            rx[i] = (y3 - y1) * (rx2[i] - x1) + (-x3 + x1) * (ry2[i] - y1)
            ry[i] = (-y2 + y1) * (rx2[i] - x1) + (x2 - x1) * (ry2[i] - y1)

        rx = 1./detJ * rx
        ry = 1./detJ * ry

        #local functions
        phi1 = 1-rx-ry
        phi2 = rx
        phi3 = ry

        #print('rx = '+ str(rx))
        #print('ry = '+ str(ry))

        #add weights
        phi = u_h[i1]*phi1 + u_h[i2]*phi2 + u_h[i3]*phi3
        #print(phi)
        r2 = tri.UniformTriRefiner(refi_t)

        t2, z2 = r2.refine_field(phi, subdiv=3)

        #ax.tricontourf(t2, z2, cmap='cool', levels=levels)
        #t2, z2 = r.refine_field(zX, subdiv=2)
        ax.tricontour(t2, z2, cmap='cool', levels=levels)



        #plt.triplot(refi_t, color='g')

    
    ax.triplot(triangulation, color='0.5', lw=0.2)
    
    #ax2 = fig.add_axes([0.82, 0.1, 0.85, 0.8])
    ax2 = fig.add_axes([0.85,0.1,0.02,0.8])
    norm = matplotlib.colors.Normalize(vmin=umin2, vmax=umax2)
    cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap='cool',
                                       norm=norm,
                                       orientation='vertical')

    plt.savefig('images/tri_sol.png', dpi=300)
    pass




def fe_plot(local_to_global, nodes_to_coordinates, borders, show_local=False):
    """Plots mesh, nodes, ... """
    
    N_fe = local_to_global.shape[0]

    fig = plt.figure()

    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlim((-1,1)) #fix constants
    ax.set_ylim((-1,1))

    X = np.zeros(shape=(N_fe, 4))
    Y = np.zeros(shape=(N_fe, 4))

    for i in range(0, N_fe):
        for j in range(0, 3):
            global_n = local_to_global[i,j]
            if j == 0:
                X[i,3] = nodes_to_coordinates[global_n, 0] 
                Y[i,3] = nodes_to_coordinates[global_n, 1]
            X[i,j] = nodes_to_coordinates[global_n, 0]
            Y[i,j] = nodes_to_coordinates[global_n, 1]

    
    for i in range(0, N_fe):
        line = matplotlib.lines.Line2D(X[i],Y[i],color="0.8")
        ax.add_line(line)



    fe_props = dict(boxstyle="round", fc="w", ec="0.8", alpha=0.9)
    node_props = dict(boxstyle="round", fc="w", ec="k", alpha=0.9)
    node_props_dir = dict(boxstyle="round", fc="w", ec="r", alpha=0.9)
    node_props_local = dict(boxstyle="round", fc="w", ec="g", alpha=0.9)

    borders2 = [b[0] for b in borders]
    borders3 = list(itertools.chain(*borders2))
    
    for i in range(0, N_fe):
        #draw fe number
        x_fe = np.sum(X[i,0:3])/3
        y_fe = np.sum(Y[i,0:3])/3
        
        #ax.text(x_fe, y_fe, str(int(local_to_global[i,j])) , ha="center", 
        #        va="center", size=5, bbox=fe_props)

        for j in range(0,3):
            if local_to_global[i,j] in borders3:
                props = node_props_dir
            else:
                props = node_props
            
            if show_local:
                ax.text(X[i,j], Y[i,j], str(j) , ha="center", 
                    va="center", size=5, bbox=node_props_local)
            else:
                ax.text(X[i,j], Y[i,j], str(int(local_to_global[i,j])) , ha="center", 
                    va="center", size=5, bbox=props)
                ax.text(X[i,j], Y[i,j], str(int(local_to_global[i,j])) , ha="center", 
                    va="center", size=5, bbox=props)


    plt.savefig('images/tri_numbers.png', dpi=300)


if __name__ == '__main__':
    
    #boundarys for pacman
    a = (4, 28, 27, 26, 25, 24, 3)
    b = (3, 23, 22, 21, 20 ,19 ,18, 2)
    c = (2, 17, 16, 15, 14, 13, 1)
    d = (1, 12, 11, 10, 9, 8, 7, 6, 0)
    e = (0, 35, 34, 33, 32, 5)
    f = (5, 31, 30, 29, 4)

    borders = [
                #(a,0),
               (b,-1),
               #(c,0),
               (d,1),
               #(e,0),
               #(f,0),
               ]

    # #boundarys for square
    # a = (3, 29, 28, 27, 26, 25, 24, 23, 22, 21, 2)
    # b = (2, 20, 19, 18, 17, 16, 15, 14, 13, 1)
    # c = (1, 12, 11, 10, 9, 8, 7, 6, 5, 4, 0)
    # d = (0, 37, 36, 35, 34, 33, 32, 31, 30, 3)

    # borders = [
    #     (c,0),
    #     (b,0),
    #     (a,0),
    # ]

    rho = lambda x: 0#math.exp(-1/(2*0.01)*((x[0]+0.4)**2+(x[1]-0.2)**2))-math.exp(-1/(2*0.01)*((x[0]+0.25)**2+(x[1]+0.35)**2))+2*math.exp(-1/(2*0.01)*((x[0]-0.2)**2+(x[1])**2))
    eps = lambda x: 1#1+(x[0]**2+x[1]**2)*1000

    local_to_global = create_local_to_global_tri('pacman_elements.txt')
    nodes_to_coordinates = create_nodes_to_coordinates_tri('pacman_coordinates.txt')
    edge_to_global, edge_to_number, edge_to_local =create_edges('pacman_edges.txt', local_to_global)

    # local_to_global = create_local_to_global_tri('square/square_elements.txt')
    # nodes_to_coordinates = create_nodes_to_coordinates_tri('square/square_coordinates.txt')
    # edge_to_global, edge_to_number, edge_to_local =create_edges('square/square_edges.txt', local_to_global)

    #nfunc = lambda x,y: 1
    func_dict = {4:lambda x,y: 0}

    neumann = {"e2g": edge_to_global, "e2n": edge_to_number, "e2l": edge_to_local, "funcs":func_dict} 


    K_h, f_h = assemble_equations(local_to_global, nodes_to_coordinates, borders, eps, rho, neumann)

    u_h = np.linalg.solve(K_h, f_h)
    print(u_h)

    contour_plot(local_to_global, nodes_to_coordinates, u_h)
    fe_plot(local_to_global, nodes_to_coordinates, borders, show_local=False)




    




