#!/usr/local/bin/python3

import os
import numpy as np

import matplotlib

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
matplotlib.rcParams['text.usetex']=True
matplotlib.rcParams['text.latex.unicode']=True

import matplotlib.pyplot as plt




def set_area(zeta, gamma, l2g, n2c, N, M, x_base, y_base, cur_node, cur_element):

    n = len(N)
    m = len(M)

    #gamma unused...

    for i in range(0,n):
        for j in range(0,m):


            x1 = x_base + sum(N[:i])
            x2 = x_base + sum(N[:i+1])
            
            y1 = y_base + sum(M[:j])
            y2 = y_base + sum(M[:j+1])

            t_l2g = [0, 0, 0, 0]


            if zeta[x1,y1] == -1:
                n2c.append([x1,y1])
                zeta[x1,y1] = cur_node
                cur_node += 1
            t_l2g[0] = zeta[x1,y1]
            

            if zeta[x2,y1] == -1:
                n2c.append([x2,y1])
                zeta[x2,y1] = cur_node
                cur_node += 1
            t_l2g[1] = zeta[x2,y1]


            if zeta[x2,y2] == -1:
                n2c.append([x2,y2])
                zeta[x2,y2] = cur_node
                cur_node += 1
            t_l2g[2] = zeta[x2,y2]


            if zeta[x1,y2] == -1:
                n2c.append([x1,y2])
                zeta[x1,y2] = cur_node
                cur_node += 1
            t_l2g[3] = zeta[x1,y2]


            l2g.append(t_l2g)
            cur_element += 1


    return cur_node, cur_element




#N = [1, 1, 2, 2, 3, 3, 4, 4, 5]
#M = [1, 1, 2, 2, 2, 3, 3, 3, 4, 4]

#upper left

xdim = 20
ydim = 20

N1 = [1 for i in range(0,xdim)]
M1 = [1 for i in range(0,ydim)]






zeta = np.zeros(shape=(sum(N1)+1,sum(M1)+1)) -1 
gamma = np.zeros(shape=zeta.shape) - 1


l2g = []
n2c = []


cur_node = 0
cur_element = 0

cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N1, M1, 0, 0, cur_node, cur_element)




nodes_to_coordinates = np.array(n2c)
local_to_global = np.array(l2g)

prop_global = dict(boxstyle="round", fc="w", ec="0.8", alpha=0.9)
prop_local = dict(boxstyle="round", fc="b", ec="0.8", alpha=0.3)


# x_len = sum(N19)+sum(N20)+ sum(N10)+sum(N9) 
# c_x_len = sum(N20)+ sum(N10)
# x_alpha = c_x_len / (1.0*x_len)

# y_len = sum(M14)+sum(M13)+sum(M12)+sum(M17)+sum(M18)+sum(M19)
# c_y_len = sum(M1)+sum(M6)
# y_alpha = c_y_len / (1.0*y_len)

# c_x_wanted = 0.1
# c_y_wanted = 0.02

# print("xlen = {}".format(x_len))
# print("c_xlen = {}".format(c_x_len))

# print("ylen = {}".format(y_len))
# print("c_ylen = {}".format(c_y_len))

print(nodes_to_coordinates[:,0].max()) #xmax 
print(nodes_to_coordinates[:,1].max()) #ymax

#print nodes_to_coordinates
xmax = nodes_to_coordinates[:,0].max() * 1.0
ymax = nodes_to_coordinates[:,1].max() * 1.0

nodes_to_coordinates = nodes_to_coordinates / 20.

print("xmax = {}".format(xmax))
print("ymax = {}".format(ymax))

#nodes_to_coordinates[:,0] = nodes_to_coordinates[:,0] * (1.0/xmax)
#scale_x = (c_x_wanted/x_alpha/xmax)
#print("scale_x = {}".format(scale_x))

#nodes_to_coordinates = nodes_to_coordinates.astype(float)
#nodes_to_coordinates[:,0] = nodes_to_coordinates[:,0] *scale_x

#new_x_max = c_x_wanted/x_alpha
#new_y_max = c_y_wanted/y_alpha

#print nodes_to_coordinates

#nodes_to_coordinates[:,1] = nodes_to_coordinates[:,1] / (nodes_to_coordinates[:,1].max()*1.0)
#nodes_to_coordinates[:,1] = nodes_to_coordinates[:,1] *(c_y_wanted/y_alpha/ymax)

#print nodes_to_coordinates
n2c = nodes_to_coordinates


N,M = zeta.shape

#boundary conditions are set here

t_edges = []
for i in range(4):
    t_edges.append([])

for i in range(N):
    for j in range(M):
        node = zeta[i,j]
        if node != -1:
            
            #outer boundary
            if j == 0:
                t_edges[0].append(node)

            if j == M-1:
                t_edges[2].append(node)

            if i == N-1: 
                t_edges[1].append(node)

            if i == 0: 
                t_edges[3].append(node)




edges = []

for e in range(0, len(t_edges)):
    if len(t_edges[e]) > 0:
        edge = np.zeros(shape=(len(t_edges[e])-1,5))

        for i in range(0,edge.shape[0]):
            if e == 0 or e == 1:
                edge[i, 0] = t_edges[e][i]
                edge[i, 1] = t_edges[e][i+1]
            else:
                edge[i, 0] = t_edges[e][i+1]
                edge[i, 1] = t_edges[e][i]

            edge[i,2] = e
            if e == 3:
                edge[i, 3] = 3
                edge[i, 4] = 0
            else:
                edge[i, 3] = e
                edge[i, 4] = e+1
    else:
        edge = np.zeros(shape=(1,2))

    edges.append(edge)

print('-1-')
print(edges[1])

print('-3-')
print(edges[3])



edge2global = np.vstack((edges[0],edges[1],edges[2],edges[3]))

print(edge2global)

draw_local = False
draw_global = False

print('nodes = ' + str(n2c.shape[0]) )
print('elements = ' + str(local_to_global.shape[0]) )





#nodes_to_coordinates = np.array(n2c)
#local_to_global = np.array(l2g)
# edges

N_nodes = nodes_to_coordinates.shape[0]

dir_set = np.zeros(shape=(1,N_nodes))
dir_value = np.zeros(shape=(1,N_nodes))

outer_bound = []
value = 0


for i in range(0,4):
    for node in list(edges[i][:,0]):
        if i == 1 or i == 3:
            dir_set[0,node] = 1
            dir_value[0,node] = value

    if i == 1 or i == 3:
        last_node = edges[i][-1,1]
        dir_set[0,last_node] = 1
        dir_value[0,last_node] = value

dir_set[0,0] = 1
dir_value[0,0] = 0


print('SUM: ' + str(np.sum(dir_set)))

boundary_table = np.vstack((dir_set, dir_value))
#print(boundary_table)
#print(dir_value)
#print(dir_set)

#print(edges[2])
#print(edges[2][::-1])

#  ------ PLOT CODE ------
draw = True

#print n2c

if draw:
    fig = plt.figure()



    for elem in range(local_to_global.shape[0]):

        l1 = local_to_global[elem, 0]
        l2 = local_to_global[elem, 1]
        l3 = local_to_global[elem, 2]
        l4 = local_to_global[elem, 3]
        if draw_local:
            plt.text(n2c[l1, 0]+0.2 , n2c[l1, 1]+0.2, "1",
                ha="center",  va="center", size=7, bbox=prop_local)

            plt.text(n2c[l2, 0]-0.2 , n2c[l2, 1]+0.2, "2",
                ha="center",  va="center", size=7, bbox=prop_local)

            plt.text(n2c[l3, 0]-0.2 , n2c[l3, 1]-0.2, "3",
                ha="center",  va="center", size=7, bbox=prop_local)

            plt.text(n2c[l4, 0]+0.2 , n2c[l4, 1]-0.2, "4",
                ha="center",  va="center", size=7, bbox=prop_local)
                
        #plt.text((n2c[l1, 0]+n2c[l2, 0])*0.5 , (n2c[l2, 1]+n2c[l3, 1])*0.5, '$'+str(elem)+'$',
        #    ha="center",  va="center", size=10)# bbox=fe_label_prop)

        x = [ n2c[l1, 0], n2c[l2, 0], n2c[l3, 0], n2c[l4, 0], n2c[l1, 0] ]
        y = [ n2c[l1, 1], n2c[l2, 1], n2c[l3, 1], n2c[l4, 1], n2c[l1, 1] ]

        plt.plot(x, y,linewidth=0.6, color='0.7')

    for i in range(0,dir_value.shape[1]):
        #if dir_value[0,i] == -1.0:
        #    plt.plot(n2c[i,0], n2c[i,1], 'ro', markersize=2)
        #    minus += 1
        if dir_set[0,i] == 1:
            #plus += 1
            plt.plot(n2c[i,0], n2c[i,1], 'ro', markersize=3)



    for e in edge2global:

        ex1 = n2c[e[0],0]       
        ey1 = n2c[e[0],1]  

        ex2 = n2c[e[1],0]       
        ey2 = n2c[e[1],1] 

        plt.arrow(ex1,ey1, ex2-ex1,ey2-ey1, head_width=0.05, head_length=0.1,length_includes_head=True)
        
        #plt.text(ex1 , ey1, str(int(e[2])), ha="center",  va="center", size=7)
    

    xlim = np.array([-0.1, 1.1])
    ylim = np.array([-0.1,1.1])
    plt.xlabel('$\mathbf{x}_1$')    
    plt.ylabel('$\mathbf{x}_2$')

    plt.gca().set_xlim(xlim)
    plt.gca().set_ylim(ylim)
    #plt.axis('equal')
    plt.savefig('images/neumann1.png', dpi=300, bbox_inches='tight')

#  ------ PLOT CODE ------

#  ------ PLOT CODE ------

# plus = 0
# minus = 0


# #plt.gca().set_xlim([-10, 100])
# #plt.gca().set_ylim([-10, 100])

# print('minus: {0} / plus: {1}'.format(minus,plus))

# plt.savefig('images/platten_mesh.pdf', dpi=300)

#  ------ PLOT CODE ------

#print nodes_to_coordinates
#print matplotlib.colors.colorConverter.to_rgba('r')


def get_nodes_to_coordinates():
    return nodes_to_coordinates

def get_edge_to_global():
    return edge2global.astype(int)

def get_local_to_global():
    return local_to_global.astype(int)

def get_boundary_table():
    return boundary_table

def get_number_of_elements():
    return local_to_global.shape[0]

def get_nodes_per_element():
    return 4

def get_number_of_nodes():
    return nodes_to_coordinates.shape[0]











