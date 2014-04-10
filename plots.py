import os

import numpy as np
import matplotlib
import matplotlib.colors as mcolors

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['text.latex.unicode']=True

import matplotlib.pyplot as plt

import helper
import fem



def generate_plot_data(info, tables, plots, u_h):
    """
    generates the x, y, z data for the plots
    """

    x = np.arange(0, plots.x_max, plots.step)
    y = np.arange(0, plots.x_max, plots.step)
    xx,yy = np.meshgrid(x, y)
    zz = np.zeros(xx.shape)

    m,n = zz.shape

    for r in range(0, info.number_of_elements):
        i1 = tables.local_to_global[r, 0]-1
        i2 = tables.local_to_global[r, 1]-1
        i3 = tables.local_to_global[r, 2]-1
        i4 = tables.local_to_global[r, 3]-1
        x1 = tables.nodes_to_coordinates[i1, 0] 
        x2 = tables.nodes_to_coordinates[i2, 0] 
        x4 = tables.nodes_to_coordinates[i4, 0] 
        y1 = tables.nodes_to_coordinates[i1, 1] 
        y2 = tables.nodes_to_coordinates[i2, 1] 
        y4 = tables.nodes_to_coordinates[i4, 1] 

        support = np.logical_and(xx >= x1, yy >= y1)

        if y4 == info.number_of_elements * 2:
            support = np.logical_and(support, yy <= y4)
        else:
            support = np.logical_and(support, yy < y4)

        if x2 == info.number_of_elements * 2:
            support = np.logical_and(support, xx <= x2)
        else:
            support = np.logical_and(support, xx < x2)

        xx2 = xx.copy()
        yy2 = yy.copy()

        _global_to_ref.elem = -1

        for i in range(0,m):
            for j in range(0,n):
                if support[i,j]:
                    xx2[i,j], yy2[i,j] = _global_to_ref(xx2[i,j],yy2[i,j],x1,x2,x4,y1,y2,y4,r)


        #local functions
        phi1 = (1-xx2) * (1-yy2)
        phi2 = xx2 * (1-yy2)
        phi3 = xx2 * yy2
        phi4 = (1-xx2) * yy2


        #add weights
        phi = u_h[i1]*phi1 + u_h[i2]*phi2 + u_h[i3]*phi3 + u_h[i4]*phi4
        phi[np.invert(support)] = 0 

        zz = zz + phi




    return xx, yy, zz

def contour_plot_ladung(info, tables, plots, u_h):
    """
    plots 
    """

    # x = np.arange(0, plots.x_max, plots.step)
    # y = np.arange(0, plots.x_max, plots.step)
    # xx,yy = np.meshgrid(x, y)
    # zz = np.zeros(xx.shape)


    plt.figure()

    for r in range(0, info.number_of_elements):



        i1 = tables.local_to_global[r, 0]
        i2 = tables.local_to_global[r, 1]
        i3 = tables.local_to_global[r, 2]
        i4 = tables.local_to_global[r, 3]
        x1 = tables.nodes_to_coordinates[i1, 0] 
        x2 = tables.nodes_to_coordinates[i2, 0] 
        x4 = tables.nodes_to_coordinates[i4, 0] 
        y1 = tables.nodes_to_coordinates[i1, 1] 
        y2 = tables.nodes_to_coordinates[i2, 1] 
        y4 = tables.nodes_to_coordinates[i4, 1] 

        x = np.arange(0, 1.2, 0.2)#plots.step)
        y = np.arange(0, 1.2, 0.2)#plots.step)
        xx,yy = np.meshgrid(x, y)



        #local functions
        phi1 = (1-xx) * (1-yy)
        phi2 = xx * (1-yy)
        phi3 = xx * yy
        phi4 = (1-xx) * yy


        #add weights
        phi = u_h[i1]*phi1 + u_h[i2]*phi2 + u_h[i3]*phi3 + u_h[i4]*phi4
        xx2 = (x2-x1)*xx + x1
        yy2 = (y4-y1)*yy + y1

        #print(phi.max())

        #levels = np.arange(u_h.min(), u_h.max()+0.01, (u_h.max()-u_h.min())/10.)
        levels = np.arange(-5.0, 5.5, 0.5)

        #maybe outside of loop
        c = mcolors.ColorConverter().to_rgb
        rvb = helper.make_colormap(
            [(1.0, 0.6, 0.4), (0.5, 0.5, 0.5), 0.5, (0.5, 0.5, 0.5), (0.0, 0.4, 0.6)])

        plt.contour(xx2, yy2, phi, cmap = plt.get_cmap('cool'), levels = levels)# , linewith=0.2) plt.get_cmap('cool')
    
    plt.xlabel('$\mathbf{x}_1$')
    plt.ylabel('$\mathbf{x}_2$')
    plt.colorbar()


    plt.savefig('images/contour_ladung.pdf', dpi=500, bbox_inches='tight')

def contour_plot_new(info, tables, plots, u_h):
    """
    plots 
    """

    # x = np.arange(0, plots.x_max, plots.step)
    # y = np.arange(0, plots.x_max, plots.step)
    # xx,yy = np.meshgrid(x, y)
    # zz = np.zeros(xx.shape)


    plt.figure()

    for r in range(0, info.number_of_elements):



        i1 = tables.local_to_global[r, 0]
        i2 = tables.local_to_global[r, 1]
        i3 = tables.local_to_global[r, 2]
        i4 = tables.local_to_global[r, 3]
        x1 = tables.nodes_to_coordinates[i1, 0] 
        x2 = tables.nodes_to_coordinates[i2, 0] 
        x4 = tables.nodes_to_coordinates[i4, 0] 
        y1 = tables.nodes_to_coordinates[i1, 1] 
        y2 = tables.nodes_to_coordinates[i2, 1] 
        y4 = tables.nodes_to_coordinates[i4, 1] 

        x = np.arange(0, 1.2, 0.2)#plots.step)
        y = np.arange(0, 1.2, 0.2)#plots.step)
        xx,yy = np.meshgrid(x, y)



        #local functions
        phi1 = (1-xx) * (1-yy)
        phi2 = xx * (1-yy)
        phi3 = xx * yy
        phi4 = (1-xx) * yy


        #add weights
        phi = u_h[i1]*phi1 + u_h[i2]*phi2 + u_h[i3]*phi3 + u_h[i4]*phi4
        xx2 = (x2-x1)*xx + x1
        yy2 = (y4-y1)*yy + y1

        #print(phi.max())

        levels = np.arange(-1.0, 1.1, 0.1)

        #maybe outside of loop
        c = mcolors.ColorConverter().to_rgb
        rvb = helper.make_colormap(
            [(0.5, 0.9, 1.0), (0.2, 0.6, 0.8), 0.5, (0.2, 0.6, 0.8), (0.0, 0.4, 0.6)])

        plt.contour(xx2, yy2, phi, cmap = plt.get_cmap('cool'), levels = levels)# , linewith=0.2) plt.get_cmap('cool')
    
    plt.xlabel('$\mathbf{x}_1$')
    plt.ylabel('$\mathbf{x}_2$')

    #plt.xlim((fem.wx1,fem.wx2))
    #plt.ylim((fem.wy1,fem.wy2))


    plt.colorbar()


    plt.savefig('images/contour_line.pdf', dpi=500, bbox_inches='tight')

def contour_plot_neumann(info, tables, plots, u_h):
    """
    plots 
    """

    # x = np.arange(0, plots.x_max, plots.step)
    # y = np.arange(0, plots.x_max, plots.step)
    # xx,yy = np.meshgrid(x, y)
    # zz = np.zeros(xx.shape)


    plt.figure()

    for r in range(0, info.number_of_elements):



        i1 = tables.local_to_global[r, 0]
        i2 = tables.local_to_global[r, 1]
        i3 = tables.local_to_global[r, 2]
        i4 = tables.local_to_global[r, 3]
        x1 = tables.nodes_to_coordinates[i1, 0] 
        x2 = tables.nodes_to_coordinates[i2, 0] 
        x4 = tables.nodes_to_coordinates[i4, 0] 
        y1 = tables.nodes_to_coordinates[i1, 1] 
        y2 = tables.nodes_to_coordinates[i2, 1] 
        y4 = tables.nodes_to_coordinates[i4, 1] 

        x = np.arange(0, 1.2, 0.2)#plots.step)
        y = np.arange(0, 1.2, 0.2)#plots.step)
        xx,yy = np.meshgrid(x, y)



        #local functions
        phi1 = (1-xx) * (1-yy)
        phi2 = xx * (1-yy)
        phi3 = xx * yy
        phi4 = (1-xx) * yy


        #add weights
        phi = u_h[i1]*phi1 + u_h[i2]*phi2 + u_h[i3]*phi3 + u_h[i4]*phi4
        xx2 = (x2-x1)*xx + x1
        yy2 = (y4-y1)*yy + y1

        #print(phi.max())

        #levels = np.arange(u_h.min(), u_h.max()+0.01, (u_h.max()-u_h.min())/10.)
        levels = np.arange(-0.25, 0.25, 0.025)

        #maybe outside of loop
        c = mcolors.ColorConverter().to_rgb
        rvb = helper.make_colormap(
            [(0.5, 0.9, 1.0), (0.2, 0.6, 0.8), 0.5, (0.2, 0.6, 0.8), (0.0, 0.4, 0.6)])

        plt.contour(xx2, yy2, phi, cmap = plt.get_cmap('cool'), levels = levels)# , linewith=0.2) plt.get_cmap('cool')
    
    plt.xlabel('$\mathbf{x}_1$')
    plt.ylabel('$\mathbf{x}_2$')

    #plt.xlim((fem.wx1,fem.wx2))
    #plt.ylim((fem.wy1,fem.wy2))


    plt.colorbar()


    plt.savefig('images/neumann.pdf', dpi=500, bbox_inches='tight')


def contour_plot(info, tables, plots, u_h):
    """Displays the data as a contour plot."""

    plt.figure()
    #plt.contourf(plots.xx, plots.yy, plots.zz)
    plt.contour(plots.xx, plots.yy, plots.zz, 15, linewidths = 0.5, colors = 'k')
    plt.pcolormesh(plots.xx, plots.yy, plots.zz, cmap = plt.get_cmap('cool'))
    plt.colorbar()


    #draw grid
    N = info.elements_per_line
    for i in range(0,N+1):
        plt.plot([i*2, i*2], [0, N*2], color='#DDDDDD', linestyle='-', linewidth=0.4)
    for i in range(0,N+1):
        plt.plot([0, N*2], [i*2, i*2], color='#DDDDDD', linestyle='-', linewidth=0.4)

    plt.savefig('images/contour.png', dpi=300)
    #plt.show()


def fe_plot(info, tables, plots):
    """Plots mesh, nodes, ... """
    
    N_fe = info.number_of_elements

    fig = plt.figure()
    #ax = plt.gca()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlim((-10, 90))
    ax.set_ylim((-10, 80))

    X = np.zeros(shape=(N_fe, 5))
    Y = np.zeros(shape=(N_fe, 5))

    for i in range(0, N_fe):


        for j in range(0, 4):
            global_n = tables.local_to_global[i,j]
            if j == 0:
                X[i,4] = tables.nodes_to_coordinates[global_n-1, 0]
                Y[i,4] = tables.nodes_to_coordinates[global_n-1, 1]
            X[i,j] = tables.nodes_to_coordinates[global_n-1, 0]
            Y[i,j] = tables.nodes_to_coordinates[global_n-1, 1]

    
    for i in range(0, N_fe):
        line = matplotlib.lines.Line2D(X[i],Y[i],color="0.8")
        ax.add_line(line)



    fe_props = dict(boxstyle="round", fc="w", ec="0.8", alpha=0.9)
    node_props = dict(boxstyle="round", fc="w", ec="k", alpha=0.9)
    node_props_dir = dict(boxstyle="round", fc="w", ec="r", alpha=0.9)

    
    for i in range(0, N_fe):
        #draw fe number
        x_fe = np.sum(X[i,0:4])/4
        y_fe = np.sum(Y[i,0:4])/4
        
        ax.text(x_fe, y_fe, str(i), ha="center", va="center", size=8, bbox=fe_props)

        for j in range(0,4):
            
            if tables.boundary_table[0,tables.local_to_global[i,j]-1] == 1:
                props = node_props_dir
            else:
                props = node_props
            ax.text(X[i,j], Y[i,j], str(int(tables.local_to_global[i,j])) , ha="center", 
                va="center", size=8, bbox=props)


    plt.savefig('images/numbers.png', dpi=300)


def _global_to_ref(x,y,x1,x2,x4,y1,y2,y4,r):
    """ Transforms global coordinates to reference coordinates"""
    J_r = np.array([[x2 - x1, x4 - x1], [y2 - y1, y4 - y1]])
    dJ_r = np.linalg.det(J_r);
    #tmp = 1/dJ_r * np.dot(np.array([[y4 - y2, -(x4-x1)], [-(y2-y1), x2-x1]]),
    #             np.array([x-x1, y-y1]))


    #optimisation
    if r > _global_to_ref.elem:
        _global_to_ref.elem += 1
    
        _global_to_ref.A = 1/dJ_r * np.array([[y4 - y2, -(x4-x1)], [-(y2-y1), x2-x1]])
        _global_to_ref.b = np.dot(_global_to_ref.A,np.array([x1, y1]))

    tmp = np.dot(_global_to_ref.A,np.array([x, y])) - _global_to_ref.b

    return tmp[0], tmp[1]


#execute run function when executed
if __name__ == '__main__':
    import fem
    fem.run()