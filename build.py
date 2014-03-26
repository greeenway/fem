#!/usr/local/bin/python3

import os
import numpy as np

import helper
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
N1 = [7, 4, 3, 2, 1, 1]
M1 = [1, 1]

N2 = [1, 1, 2, 2, 3, 3, 4, 4, 5]
M2 = M1[:]

N3 = N2[:]
M3 = [1, 1]

N4 = N2[:]
M4 = [1, 1, 2, 2, 2, 3, 3, 4, 4, 5]

N5 = N1[:]
M5 = M4[:]

#upper right
N11 = N1[::-1] #reversed
M11 = M1[:]

N12 = N2[::-1]
M12 = M1[:]

N13 = N12[:]
M13 = M3[:]

N14 = N12[:]
M14 = M4[:]

N15 = N5[::-1]
M15 = M5[:]

#lower right
N6 = N1
M6 = M1[::-1]

N7 = N2
M7 = M2[::-1]

N8 = N7
M8 = M3[::-1]

N9 = N7
M9 = M4[::-1]

N10 = N6
M10 = M9

#lower left
N16 = N11
M16 = M11[::-1]

N17 = N12
M17 = M12[::-1]

N18 = N17
M18 = M13[::-1]

N19 = N17
M19 = M14[::-1]

N20 = N10[::-1]
M20 = M19




zeta = np.zeros(shape=(sum(N1)+sum(N2)+sum(N11)+sum(N12)+1,sum(M2)+sum(M3)+sum(M4)+sum(M9)+sum(M8)+sum(M7)+1)) -1 
gamma = np.zeros(shape=zeta.shape) - 1


l2g = []
n2c = []


cur_node = 0
cur_element = 0

cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N1, M1, sum(N11)+sum(N12), sum(M9)+sum(M8)+sum(M7), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N2, M2, sum(N11)+sum(N12)+sum(N1), sum(M9)+sum(M8)+sum(M7), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N3, M3, sum(N11)+sum(N12)+sum(N1), sum(M2)+sum(M9)+sum(M8)+sum(M7), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N4, M4, sum(N11)+sum(N12)+sum(N1), sum(M2)+sum(M3)+sum(M9)+sum(M8)+sum(M7), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N5, M5, sum(N11)+sum(N12), sum(M2)+sum(M3)+sum(M9)+sum(M8)+sum(M7), cur_node, cur_element)

cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N11, M11, sum(N12), sum(M9)+sum(M8)+sum(M7), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N12, M12, 0, sum(M9)+sum(M8)+sum(M7), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N13, M13, 0, sum(M12)+sum(M9)+sum(M8)+sum(M7), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N14, M14, 0, sum(M12)+sum(M13)+sum(M9)+sum(M8)+sum(M7), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N15, M15, sum(N14), sum(M12)+sum(M13)+sum(M9)+sum(M8)+sum(M7), cur_node, cur_element)

cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N6, M6, sum(N11)+sum(N12), sum(M9)+sum(M8), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N7, M7, sum(N11)+sum(N12)+sum(N6), sum(M9)+sum(M8), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N8, M8, sum(N11)+sum(N12)+sum(N6), sum(M9), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N9, M9, sum(N11)+sum(N12)+sum(N6), 0, cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N10, M10, sum(N11)+sum(N12), 0, cur_node, cur_element)

cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N16, M16, sum(N17), sum(M19)+sum(M18), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N17, M17, 0, sum(M19)+sum(M18), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N18, M18, 0, sum(M19), cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N19, M19, 0, 0, cur_node, cur_element)
cur_node, cur_element = set_area(zeta, gamma, l2g, n2c, N20, M20, sum(N19), 0, cur_node, cur_element)



nodes_to_coordinates = np.array(n2c)
local_to_global = np.array(l2g)

prop_global = dict(boxstyle="round", fc="w", ec="0.8", alpha=0.9)
prop_local = dict(boxstyle="round", fc="b", ec="0.8", alpha=0.3)


x_len = sum(N19)+sum(N20)+ sum(N10)+sum(N9) 
c_x_len = sum(N20)+ sum(N10)
x_alpha = c_x_len / (1.0*x_len)

y_len = sum(M14)+sum(M13)+sum(M12)+sum(M17)+sum(M18)+sum(M19)
c_y_len = sum(M1)+sum(M6)
y_alpha = c_y_len / (1.0*y_len)

c_x_wanted = 0.1
c_y_wanted = 0.02

print("xlen = {}".format(x_len))
print("c_xlen = {}".format(c_x_len))

print("ylen = {}".format(y_len))
print("c_ylen = {}".format(c_y_len))

print(nodes_to_coordinates[:,0].max()) #xmax 
print(nodes_to_coordinates[:,1].max()) #ymax

#print nodes_to_coordinates
xmax = nodes_to_coordinates[:,0].max() * 1.0
ymax = nodes_to_coordinates[:,1].max() * 1.0

print("xmax = {}".format(xmax))
print("ymax = {}".format(ymax))

#nodes_to_coordinates[:,0] = nodes_to_coordinates[:,0] * (1.0/xmax)
scale_x = (c_x_wanted/x_alpha/xmax)
print("scale_x = {}".format(scale_x))



nodes_to_coordinates = nodes_to_coordinates.astype(float)

x1 = (sum(N17)*1.0)*scale_x
print('x1 = {}'.format(x1))

y1 = ((sum(M19)+sum(M18))*1.0)*(c_y_wanted/y_alpha/ymax)
print('y1 = {}'.format(y1))


x2 = ((sum(N17)+sum(N16)+sum(N6))*1.0)*scale_x
print('x2 = {}'.format(x1))

y2 = ((sum(M19)+sum(M18)+sum(M12)+sum(M17))*1.0)*(c_y_wanted/y_alpha/ymax)
print('y2 = {}'.format(y1))

important_points = np.array([[x1, y1],[x2, y2]])
helper.write_matrix_to_file(important_points, 'data/important_points.txt')

wx1 = ((sum(N17)+sum(N16)+sum(N6)-15)*1.0)*scale_x


wy1 = ((sum(M19))*1.0)*(c_y_wanted/y_alpha/ymax)


wx2 = ((sum(N17)+sum(N16)+sum(N6)+sum(N7)-15)*1.0)*scale_x

wy2 = ((sum(M19)+sum(M18)+sum(M12)+sum(M17)+sum(M12))*1.0)*(c_y_wanted/y_alpha/ymax)

important_window = np.array([[wx1, wy1],[wx2, wy2]])
helper.write_matrix_to_file(important_window, 'data/window.txt')

nodes_to_coordinates[:,0] = nodes_to_coordinates[:,0] *scale_x




new_x_max = c_x_wanted/x_alpha
new_y_max = c_y_wanted/y_alpha

#print nodes_to_coordinates

#nodes_to_coordinates[:,1] = nodes_to_coordinates[:,1] / (nodes_to_coordinates[:,1].max()*1.0)
nodes_to_coordinates[:,1] = nodes_to_coordinates[:,1] *(c_y_wanted/y_alpha/ymax)

#print nodes_to_coordinates
n2c = nodes_to_coordinates


N,M = zeta.shape

t_edges = []
for i in range(12):
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

            #lower electrode
            if j == sum(M19) and i >= sum(N18) and i <= sum(N18)+sum(N16)+sum(N6):
                t_edges[4].append(node)

            if j >= sum(M9) and j <= sum(M9)+sum(M8) and i == sum(N18)+sum(N16)+sum(N6):
                t_edges[5].append(node)

            if j == sum(M19)+sum(M18) and i >= sum(N18) and i <= sum(N18)+sum(N16)+sum(N6):
                t_edges[6].append(node)

            if j >= sum(M9) and j <= sum(M9)+sum(M8) and i == sum(N18):
                t_edges[7].append(node)

            #upper electride
            if j == sum(M19)+sum(M18)+sum(M17)+sum(M12) and i >= sum(N18) and i <= sum(N18)+sum(N16)+sum(N6):
                t_edges[8].append(node)

            if j >= sum(M19)+sum(M18)+sum(M17)+sum(M12) and j <= sum(M19)+sum(M18)+sum(M17)+sum(M12)+sum(M13) and i == sum(N18)+sum(N16)+sum(N6):
                t_edges[9].append(node)

            if j == sum(M19)+sum(M18)+sum(M17)+sum(M12)+sum(M13) and i >= sum(N18) and i <= sum(N18)+sum(N16)+sum(N6):
                t_edges[10].append(node)

            if j >= sum(M19)+sum(M18)+sum(M17)+sum(M12) and j <= sum(M19)+sum(M18)+sum(M17)+sum(M12)+sum(M13) and i == sum(N18):
                t_edges[11].append(node)

edges = []

for e in range(0, len(t_edges)):
    if len(t_edges[e]) > 0:
        edge = np.zeros(shape=(len(t_edges[e])-1,2))

        for i in range(0,edge.shape[0]):
            edge[i, 0] = t_edges[e][i]
            edge[i, 1] = t_edges[e][i+1]
    else:
        edge = np.zeros(shape=(1,2))

    edges.append(edge)


#print(lower)

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
outer_value = 0
lower_value = -1
higher_value = 1

for i in range(0,4):
    for node in list(edges[i][:,0]):
        dir_set[0,node] = 1
        dir_value[0,node] = outer_value

    last_node = edges[i][-1,1]
    dir_set[0,last_node] = 1
    dir_value[0,last_node] = outer_value


for i in range(4,8):
    for node in list(edges[i][:,0]):
        dir_set[0,node] = 1
        dir_value[0,node] = lower_value
    last_node = edges[i][-1,1]
    dir_set[0,last_node] = 1
    dir_value[0,last_node] = lower_value

for i in range(8,12):
    for node in list(edges[i][:,0]):
        dir_set[0,node] = 1
        dir_value[0,node] = higher_value
    last_node = edges[i][-1,1]
    dir_set[0,last_node] = 1
    dir_value[0,last_node] = higher_value


boundary_table = np.vstack((dir_set, dir_value))

#  ------ PLOT CODE ------
fig = plt.figure()

for elem in range(local_to_global.shape[0]):

    l1 = local_to_global[elem, 0]
    l2 = local_to_global[elem, 1]
    l3 = local_to_global[elem, 2]
    l4 = local_to_global[elem, 3]
    # if draw_local:
    #     plt.text(n2c[l1, 0]+0.2 , n2c[l1, 1]+0.2, "1",
    #         ha="center",  va="center", size=7, bbox=prop_local)

    #     plt.text(n2c[l2, 0]-0.2 , n2c[l2, 1]+0.2, "2",
    #         ha="center",  va="center", size=7, bbox=prop_local)

    #     plt.text(n2c[l3, 0]-0.2 , n2c[l3, 1]-0.2, "3",
    #         ha="center",  va="center", size=7, bbox=prop_local)

    #     plt.text(n2c[l4, 0]+0.2 , n2c[l4, 1]-0.2, "4",
    #         ha="center",  va="center", size=7, bbox=prop_local)
            
    # plt.text((n2c[l1, 0]+n2c[l2, 0])*0.5 , (n2c[l2, 1]+n2c[l3, 1])*0.5, str(elem),
    #     ha="center",  va="center", size=2)# bbox=fe_label_prop)

    x = [ n2c[l1, 0], n2c[l2, 0], n2c[l3, 0], n2c[l4, 0], n2c[l1, 0] ]
    y = [ n2c[l1, 1], n2c[l2, 1], n2c[l3, 1], n2c[l4, 1], n2c[l1, 1] ]

    plt.plot(x, y,linewidth=0.6, color='0.7')

plt.plot(x1, y1, 'ro')
plt.plot(x2, y2, 'bo')

xlim = np.array([0, new_x_max])
ylim = np.array([0, new_y_max])
plt.xlabel('$\mathbf{x}_1$')    
plt.ylabel('$\mathbf{x}_2$')

#plt.gca().set_xlim(xlim)
#plt.gca().set_ylim(ylim)
plt.xlim(wx1,wx2)
plt.ylim(wy1,wy2)
#plt.axis('equal')
plt.savefig('images/platten_mesh.pdf', dpi=500, bbox_inches='tight')

#  ------ PLOT CODE ------

#  ------ PLOT CODE ------

# plus = 0
# minus = 0

# for i in range(0,dir_value.shape[1]):
#     if dir_value[0,i] == -1.0:
#         plt.plot(n2c[i,0], n2c[i,1], 'ro', markersize=2)
#         minus += 1
#     if dir_value[0,i] == 1.0:
#         plus += 1
#         plt.plot(n2c[i,0], n2c[i,1], 'bo', markersize=2)

# #plt.gca().set_xlim([-10, 100])
# #plt.gca().set_ylim([-10, 100])

# print('minus: {0} / plus: {1}'.format(minus,plus))

# plt.savefig('images/platten_mesh.pdf', dpi=300)

#  ------ PLOT CODE ------

print nodes_to_coordinates
print matplotlib.colors.colorConverter.to_rgba('r')


def get_nodes_to_coordinates():
    return nodes_to_coordinates


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


print('sum = ' + str( sum(N12)  ))
# print('sum = ' + str( sum(M19)+sum(M18)+sum(M17)+sum(M12)))










