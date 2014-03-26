#!/usr/local/bin/python3

import numpy as np
#import matplotlib
import matplotlib.pyplot as plt

import numeric
import helper
from helper import Timer
import fem
import assembly

#def rho(x): return 0#math.exp(-0.5/4**2*((x[0]-4)**2+(x[1]-12)**2))*1e-10-math.exp(-0.5/4**2*((x[0]-20)**2+(x[1]-12)**2))*1e-10
#def eps(x): return 1#return if ((x[0] > x1 ) and (x[0] < x2) and (x[1] > y1) and (x[1] < y2))#1

with Timer('importing solution'):
    nodes_to_coordinates = helper.read_matrix_from_file('data/nodes_to_coordinates_platten.txt')
    

    local_to_global = helper.read_matrix_from_file('data/local_to_global_platten.txt')
    boundary_table = helper.read_matrix_from_file('data/boundary_table_platten.txt')

    #print nodes_to_coordinates

    u_h = helper.read_matrix_from_file('data/u_h_platten.txt')


with Timer('Calculation 1'):
    energy = 0

     # phi1 = (1-xx) * (1-yy)
     # phi2 = xx * (1-yy)
     # phi3 = xx * yy
     # phi4 = (1-xx) * yy

     #grad phi1 = [y-1, x-1]
     #grad phi2 = [1-y, -x]
     #grad phi3 = [y, x]
     #grad phi4 = [-y, 1-x ]

    N_fe = local_to_global.shape[0]

    eps = 8.854e-12
    energies = np.zeros(shape=(N_fe,1))

    for r in range(0, N_fe):
        i1 = local_to_global[r, 0]
        i2 = local_to_global[r, 1]
        i3 = local_to_global[r, 2]
        i4 = local_to_global[r, 3]


        x1 = nodes_to_coordinates[i1, 0] 
        x2 = nodes_to_coordinates[i2, 0] 
        x4 = nodes_to_coordinates[i4, 0] 
        y1 = nodes_to_coordinates[i1, 1] 
        y2 = nodes_to_coordinates[i2, 1] 
        y4 = nodes_to_coordinates[i4, 1] 

        J_r = np.array([[x2 - x1, x4 - x1], [y2 - y1, y4 - y1]])
        dJ_r = np.linalg.det(J_r)
        adJ_r = abs(dJ_r)

        #Jit = np.linalg.inv(J_r).transpose()

        u1 = u_h[i1,0]
        u2 = u_h[i2,0]
        u3 = u_h[i3,0]
        u4 = u_h[i4,0]

        #if fem.eps((0.1,0.15)) > 8.854e-12:
        #    print 'larger'

        #print fem.eps((0.1,0.15))

        #g1 = lambda x,y: ((-y + x) * y1 - y4 + y4 * y + y2 - y2 * x)**2 + ((-y + x) * x1 - x4 + x4 * y + x2 - x2 * x)**2
        #g2 = lambda x,y: (-2 * y2 * y1 - 2 * x2 * x1 + x2 * x2 + x1 * x1 + y2 * y2 + y1 * y1) * x * x - 2 * (-1 + y) * (x1 * x1 + (-x2 - x4) * x1 + y1 * y1 + (-y2 - y4) * y1 + x4 * x2 + y4 * y2) * x + (-2 * y4 * y1 - 2 * x4 * x1 + x4 * x4 + x1 * x1 + y4 * y4 + y1 * y1) * ((-1 + y)**2)

        #g3 = lambda x,y: ((-y + x) * y1 - y2 * x + y4 * y)**2 + ((-y + x) * x1 - x2 * x + x4 * y)**2
        #g4 = lambda x,y: y * y * (-2 * y4 * y1 - 2 * x4 * x1 + x4 * x4 + x1 * x1 + y4 * y4 + y1 * y1) - 2 * (-1 + x) * (x1 * x1 + (-x2 - x4) * x1 + y1 * y1 + (-y2 - y4) * y1 + x4 * x2 + y4 * y2) * y + (-2 * y2 * y1 - 2 * x2 * x1 + x2 * x2 + x1 * x1 + y2 * y2 + y1 * y1) * (-1 + x)**2

        #divide bei 1/detJ^2
        #integrand = lambda xi: eps*(u1**2*g1(xi[0],xi[1])) + u2**2*g2(xi[0],xi[1]) + u3**2*g3(xi[0],xi[1]) + u4**2*g4(xi[0],xi[1]) )

        #integrand = lambda x,y: eps*(u1**2*g1(x,y) + u2**2*g2(x,y) + u3**2*g3(x,y) + u4**2*g4(x,y) )
        
        #integrand = lambda x,y: eps* ( ((y4-y1)*(u1*(-1+y)+u2*(1-y)+u3*y-u4*y)+(y1-y2)*(u1*(x-1)-u2*x+u3*x+u4*(1-x)))**2+((x1-x4)*(u1*(y-1)+u2*(1-y)+u3*y-u4*y)+(x2-x1)*(u1*(x-1)-u2*x+u3*x+u4*(1-x))))**2   )

        u = [u1, u2, u3, u4]
        integrand = lambda x,y: fem.eps(assembly._ref_to_global(x1, x2, x4, y1, y2, y4, (x,y)))*((-y4 * u[0] + y4 * u[0] * y + y4 * u[1] - y4 * u[1] * y + y4 * u[2] * y - y4 * u[3] * y - y1 * u[0] * y - y1 * u[1] + y1 * u[1] * y - y1 * u[2] * y + y1 * u[3] * y + y2 * u[0] - y2 * u[0] * x + y2 * u[1] * x - y2 * u[2] * x - y2 * u[3] + y2 * u[3] * x + y1 * u[0] * x - y1 * u[1] * x + y1 * u[2] * x + y1 * u[3] - y1 * u[3] * x)**2 + (-x4 * u[0] + x4 * u[0] * y + x4 * u[1] - x4 * u[1] * y + x4 * u[2] * y - x4 * u[3] * y - x1 * u[0] * y - x1 * u[1] + x1 * u[1] * y - x1 * u[2] * y + x1 * u[3] * y + x2 * u[0] - x2 * u[0] * x + x2 * u[1] * x - x2 * u[2] * x - x2 * u[3] + x2 * u[3] * x + x1 * u[0] * x - x1 * u[1] * x + x1 * u[2] * x + x1 * u[3] - x1 * u[3] * x)**2)

        energy_i = 0

        if True:#x1 >= 25 and x2 <= 57 and y1 >= 27 and y4 <= 33: #only between electrodes
        #if x1 >= 25 and x2 <= 57 and y1 >= 19 and y4 <= 23: #only between electrodes
            energy_i = adJ_r/dJ_r**2 * numeric.gauss_2D_3(integrand)
        


        #print([x1, x2, y1, y4])
        #print(energy_i)

        #integrand = lambda x,y: eps*(u1**2*g1(x,y) + u2**2*g2(x,y) + u3**2*g3(x,y) + u4**2*g4(x,y) ) #wrong...:(

        energies[r,0] = energy_i

        #if x1 >= 25 and x2 <= 57 and y1 >= 19 and y4 <= 23: #only between electrodes
        #     energy += adJ_r/dJ_r**2 * numeric.gauss_2D_3(integrand)

        #print(adJ_r/dJ_r**2 * numeric.gauss_2D_3(integrand))
        energy += energy_i




    C_ = energy/(1-(-1))**2

 
    print('C_ = ' + str(C_))
    print('C_inf = ' + str(eps*10*0.1/0.02))
#  ------ PLOT CODE ------    
if False:
    with Timer('plotting'):
        

    # 
        helper.write_matrix_to_file(energies, 'data/energies.txt')

        emax = energies.max()
        emin = energies.min()

        e_tmp = energies #/ emax


        fig = plt.figure() 
        ax = plt.gca()

        from matplotlib.patches import Ellipse, Polygon

        ax.fill([-10,90,90,-10],[-10,-10,70,70], fill=True, color='b')

        for r in range(0, N_fe):
            i1 = local_to_global[r, 0]
            i2 = local_to_global[r, 1]
            i3 = local_to_global[r, 2]
            i4 = local_to_global[r, 3]
            x1 = nodes_to_coordinates[i1, 0] 
            x2 = nodes_to_coordinates[i2, 0]
            x3 = nodes_to_coordinates[i3, 0]  
            x4 = nodes_to_coordinates[i4, 0] 
            y1 = nodes_to_coordinates[i1, 1] 
            y2 = nodes_to_coordinates[i2, 1] 
            y3 = nodes_to_coordinates[i3, 1] 
            y4 = nodes_to_coordinates[i4, 1] 

            e_tmp[r][0] = e_tmp[r][0] /  ( (x2-x1)*(y4-y1) )

        emax2 = e_tmp.max()
        e_tmp = e_tmp / emax2

        for r in range(0, N_fe):
            i1 = local_to_global[r, 0]
            i2 = local_to_global[r, 1]
            i3 = local_to_global[r, 2]
            i4 = local_to_global[r, 3]
            x1 = nodes_to_coordinates[i1, 0] 
            x2 = nodes_to_coordinates[i2, 0]
            x3 = nodes_to_coordinates[i3, 0]  
            x4 = nodes_to_coordinates[i4, 0] 
            y1 = nodes_to_coordinates[i1, 1] 
            y2 = nodes_to_coordinates[i2, 1] 
            y3 = nodes_to_coordinates[i3, 1] 
            y4 = nodes_to_coordinates[i4, 1] 

            ax.fill([x1,x2,x3,x4],[y1,y2,y3,y4], fill=True, color=str(e_tmp[r][0]))

            # if x1 >= 25 and x2 <= 57 and y1 >= 27 and y4 <= 33:
            #     plt.plot(x4,y4,'or')
        plt.gca().set_xlim([-0.1, 1.1])
        plt.gca().set_ylim([-0.1, 1.1])
        plt.savefig('images/energy_distri.png', dpi=300)
#  ------ PLOT CODE ------ 









#execute run function when executed
if __name__ == '__main__':
    pass