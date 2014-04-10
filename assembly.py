import numpy as np
import math

import fem
import numeric


linear_trial_funcs = [
    lambda xi1, xi2: (1-xi1)*(1-xi2),
    lambda xi1, xi2: xi1*(1-xi2),
    lambda xi1, xi2: xi1*xi2,
    lambda xi1, xi2: (1-xi1)*xi2,
]

linear_trial_funcs_1D = [
    lambda xi: (1-xi),
    lambda xi: xi,
]

dn_funcs = [
    lambda xi: -math.sin(math.pi*xi[0]),
    lambda xi: 1,
    lambda xi: math.sin(math.pi*xi[0]),
    lambda xi: 1,
]



def assemble_equations(info, tables):
    """
    Assembles the final system of linear equations to be solved.
    Outputs a Matrix K_h and the right hand side vector f_h
    """

    N_nodes = info.number_of_nodes
    N_fe = info.number_of_elements

    K_h = np.zeros(shape=(N_nodes, N_nodes))
    f_h = np.zeros(shape=(N_nodes,1))


    for r in range(0,N_fe):
        #print(N_fe)
        i1 = tables.local_to_global[r, 0]
        i2 = tables.local_to_global[r, 1]
        i4 = tables.local_to_global[r, 3]
        x1 = tables.nodes_to_coordinates[i1, 0] 
        x2 = tables.nodes_to_coordinates[i2, 0] 
        x4 = tables.nodes_to_coordinates[i4, 0] 
        y1 = tables.nodes_to_coordinates[i1, 1] 
        y2 = tables.nodes_to_coordinates[i2, 1] 
        y4 = tables.nodes_to_coordinates[i4, 1] 

        K_r = _calculate_K_r(x1, x2, x4, y1, y2, y4, fem.eps) #TESTING #np.zeros(shape=(4,4))
        f_r = _calculate_f_r(x1, x2, x4, y1, y2, y4, fem.rho) #np.zeros(shape=(4,1))#

        for alpha in range(0,4):
            i = tables.local_to_global[r,alpha]# - 1
            
            f_h[i,0] = f_h[i,0] + f_r[alpha,0]
            
            for beta in range(0,4):
                j = tables.local_to_global[r,beta]# - 1

                K_h[i,j] = K_h[i,j] + K_r[alpha,beta] 

    # --- neumann ---
    if False:
        N_edges = tables.edge_to_global.shape[0]
        e_tab = tables.edge_to_global

        #print N_edges
        for e in range(0,N_edges):
            #calculate f^(e)

            f_e = _calculate_f_e(info, tables, e )

            for alpha in range(0,2):

                f_h[ e_tab[e, alpha] ] = f_e[alpha]



    # --- neumann ---


    #dirichlet boundary
    b_tmp = tables.boundary_table[0,:]*tables.boundary_table[1,:]

    for i in range(0, N_nodes):
        for j in range(0, N_nodes):
            if tables.boundary_table[0,i] == 0:
                f_h[i] = f_h[i] - K_h[i,j] * b_tmp[j]

    for i in range(0,N_nodes):
        for j in range(0,N_nodes):
            if tables.boundary_table[0,i] == 1 or \
                tables.boundary_table[0,j] == 1:
                if i == j:
                    K_h[i,j] = 1
                else:
                    K_h[i,j] = 0

        if tables.boundary_table[0,i] == 1:
            f_h[i] = tables.boundary_table[1,i]

    equ = fem.Equations()
    equ.K_h = K_h
    equ.f_h = f_h

    return equ

def _calculate_f_e(info, tables , e):
    f_e = np.zeros(shape=(2,1))

    e_tab = tables.edge_to_global
    f1 = e_tab[e, 3] #numbers of functions in T^(r) -> linear_trial_funcs[f1]
    f2 = e_tab[e, 4]

    i1 = e_tab[e, 0]
    i2 = e_tab[e, 1]
    nr = e_tab[e, 2]

    x1 = tables.nodes_to_coordinates[i1, 0]
    y1 = tables.nodes_to_coordinates[i1, 1]

    x2 = tables.nodes_to_coordinates[i2, 0]
    y2 = tables.nodes_to_coordinates[i2, 1]

    length = math.sqrt((x2-x1)**2+(y2-y1)**2)

    integ2 = lambda xi: dn_funcs[nr](_ref_e_to_global(x1, x2, y1, y2, xi))*linear_trial_funcs_1D[1](xi) * length
    integ1 = lambda xi: dn_funcs[nr](_ref_e_to_global(x1, x2, y1, y2, xi))*linear_trial_funcs_1D[0](xi) * length

    f_e[0] = numeric.gauss_1D_6(integ1)
    f_e[1] = numeric.gauss_1D_6(integ2)


    return f_e

def _ref_e_to_global(x1, x2, y1, y2, xi):

    tx1 = (x2-x1)*xi + x1
    tx2 = (y2-y1)*xi + y1

    return (tx1,tx2)

def _calculate_K_r(x1, x2, x4, y1, y2, y4, eps): 
    """
    _calculate_K_r
    """
    J_r = np.array([[x2 - x1, x4 - x1], [y2 - y1, y4 - y1]])
    dJ_r = np.linalg.det(J_r);

    J11 =  (y4 - y1)**2 + (-x4 + x1)**2
    J12 =(y4 - y1) * (-y2 + y1) + (-x4 + x1) * (x2 - x1)
    J21 =(y4 - y1) * (-y2 + y1) + (-x4 + x1) * (x2 - x1)
    J22 =  (-y2 + y1)**2 + (x2 - x1)**2;


    # K_r = np.zeros(shape=(4,4)) #remove constants?
    # K_r[0,0] = 1/3. * J11 + 1/4. * J12 + 1/4. *J21 + 1/3.*J22
    # K_r[0,1] = -1/3. * J11 + 1/4. * J12 - 1/4. * J21 + 1/6.*J22
    # K_r[0,2] = -1/6. * J11 - 1/4. * J12 - 1/4. * J21 - 1/6.* J22
    # K_r[0,3] = 1/6. *J11 - 1/4.*J12 + 1/4. *J21 - 1/3. *J22

    # K_r[1,0] = -1/3.*J11 - 1/4.*J12 + 1/4.*J21 + 1/6.*J22
    # K_r[1,1] = 1/3.*J11 - 1/4.*J12 - 1/4.*J21 + 1/3.*J22
    # K_r[1,2] = 1/6. *J11 + 1/4.*J12 - 1/4.*J21 - 1/3.*J22
    # K_r[1,3] = -1/6.*J11 + 1/4.*J12 + 1/4.*J21 - 1/6.*J22

    # K_r[2,0] = -1/6.*J11-1/4.*J12 -1/4.*J21 - 1/6.*J22
    # K_r[2,1] = 1/6.*J11-1/4.*J12+1/4.*J21 - 1/3.*J22
    # K_r[2,2] = 1/3.*J11 + 1/4.*J12 + 1/4.*J21 + 1/3.*J22
    # K_r[2,3] = -1/3.*J11 + 1/4.*J12 - 1/4.*J21 + 1/6.*J22

    # K_r[3,0] = 1/6.*J11 + 1/4.*J12-1/4.*J21 - 1/3.*J22
    # K_r[3,1] = -1/6.*J11 + 1/4.*J12+1/4.*J21 - 1/6.*J22
    # K_r[3,2] = -1/3.*J11 - 1/4.*J12 + 1/4.*J21 + 1/6.*J22
    # K_r[3,3] = 1/3.*J11 -1/4.*J12 - 1/4.*J21 + 1/3.*J22

    # K_r = 1 * abs(dJ_r)/ (dJ_r**2) *  K_r;

    K_r2 = np.zeros(shape=(4,4)) #for numeric
    #K_r2[0,0] = numeric.gauss_2D_3(lambda x,y: ((J11 * (-1 + y) + J12 * (-1 + x)) * (-1 + y) + (J21 * (-1 + y) + J22 * (-1 + x)) * (-1 + x)))


    #eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y)))

    K_r2[0,0] = numeric.gauss_2D_3(lambda x,y: ((J11 * (-1 + y) + J12 * (-1 + x)) * (-1 + y) + (J21 * (-1 + y) + J22 * (-1 + x)) * (-1 + x)) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[0,1] = numeric.gauss_2D_3(lambda x,y: ((J11 * (1 - y) - J12 * x) * (-1 + y) + (J21 * (1 - y) - J22 * x) * (-1 + x)) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[0,2] = numeric.gauss_2D_3(lambda x,y: ((J11 * y + J12 * x) * (-1 + y) + (J21 * y + J22 * x) * (-1 + x)) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[0,3] = numeric.gauss_2D_3(lambda x,y: ((-J11 * y + J12 * (1 - x)) * (-1 + y) + (-J21 * y + J22 * (1 - x)) * (-1 + x)) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[1,0] = numeric.gauss_2D_3(lambda x,y: ((J11 * (-1 + y) + J12 * (-1 + x)) * (1 - y) - (J21 * (-1 + y) + J22 * (-1 + x)) * x) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[1,1] = numeric.gauss_2D_3(lambda x,y: ((J11 * (1 - y) - J12 * x) * (1 - y) - (J21 * (1 - y) - J22 * x) * x) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[1,2] = numeric.gauss_2D_3(lambda x,y: ((J11 * y + J12 * x) * (1 - y) - (J21 * y + J22 * x) * x) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[1,3] = numeric.gauss_2D_3(lambda x,y: ((-J11 * y + J12 * (1 - x)) * (1 - y) - (-J21 * y + J22 * (1 - x)) * x) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[2,0] = numeric.gauss_2D_3(lambda x,y: ((J11 * (-1 + y) + J12 * (-1 + x)) * y + (J21 * (-1 + y) + J22 * (-1 + x)) * x) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[2,1] = numeric.gauss_2D_3(lambda x,y: ((J11 * (1 - y) - J12 * x) * y + (J21 * (1 - y) - J22 * x) * x) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[2,2] = numeric.gauss_2D_3(lambda x,y: ((J11 * y + J12 * x) * y + (J21 * y + J22 * x) * x) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[2,3] = numeric.gauss_2D_3(lambda x,y: ((-J11 * y + J12 * (1 - x)) * y + (-J21 * y + J22 * (1 - x)) * x) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[3,0] = numeric.gauss_2D_3(lambda x,y: (-(J11 * (-1 + y) + J12 * (-1 + x)) * y + (J21 * (-1 + y) + J22 * (-1 + x)) * (1 - x)) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[3,1] = numeric.gauss_2D_3(lambda x,y: (-(J11 * (1 - y) - J12 * x) * y + (J21 * (1 - y) - J22 * x) * (1 - x)) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[3,2] = numeric.gauss_2D_3(lambda x,y: (-(J11 * y + J12 * x) * y + (J21 * y + J22 * x) * (1 - x)) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))
    K_r2[3,3] = numeric.gauss_2D_3(lambda x,y: (-(-J11 * y + J12 * (1 - x)) * y + (-J21 * y + J22 * (1 - x)) * (1 - x)) * eps(_ref_to_global(x1, x2, x4, y1, y2, y4, (x,y))))


    K_r2 = 1./ abs(dJ_r) *  K_r2;
    #print(K_r-K_r2)
    return K_r2


def _calculate_f_r(x1, x2, x4, y1, y2, y4, rho):
    """
    _calculate_f_r elementwise "lastvektoren"
    """
    J_r = np.array([[x2 - x1, x4 - x1], [y2 - y1, y4 - y1]])
    a_dJ_r = abs(np.linalg.det(J_r))
    
    N = 4
    f_r = np.zeros(shape=(N,1))
    phi = linear_trial_funcs
    
    for i in range(0,4):
        
        integrand = lambda xi1, xi2: fem.rho(_ref_to_global(x1, x2, x4, y1, y2, y4, (xi1,xi2)))*phi[i](xi1,xi2)*a_dJ_r
        f_r[i,0] = numeric.gauss_2D_3(integrand)

    return f_r
    

def _ref_to_global(x1, x2, x4, y1, y2, y4, xi):
    xi1 = xi[0]
    xi2 = xi[1]
    tmp = np.dot(np.array([[x2-x1, x4-x1], [y2-y1, y4-y1]]), np.array([xi1, xi2])) + np.array([x1,y1]);
    return tmp[0],tmp[1]

#execute run function when executed
if __name__ == '__main__':
    import fem
    fem.run()













