#!/usr/local/bin/python3

import math

import preparation
import assembly
import solve
import plots
import helper

from helper import Timer

ELEMENTS_PER_LINE = 12
NODES_PER_ELEMENT = 4
STEP_SIZE = 0.2

points = helper.read_matrix_from_file('data/important_points.txt')
print('x1='+str(points[0,0])+'y1='+str(points[0,1]))
x1 = points[0,0]
y1 = points[0,1]

x2 = points[1,0]
y2 = points[1,1]

window = helper.read_matrix_from_file('data/window.txt')
wx1 = window[0,0]
wy1 = window[0,1]

wx2 = window[1,0]
wy2 = window[1,1]

def rho(x): return 0#math.exp(-0.5/4**2*((x[0]-7)**2+(x[1]-10)**2))-math.exp(-0.5/4**2*((x[0]-13)**2+(x[1]-10)**2))

eps0 = 8.854*10**(-12)

#def eps(x): 
#    return eps0

def eps(x): 
    if  ( (x[0] > x1) and (x[0] < x2) and (x[1] > y1) and (x[1] < y2)):
       return eps0*10
    else: 
       return eps0 #1 #return 1#+math.exp(-0.5/4**2*((x[1]-12)**2))*100


#and (x[0] < x2) and (x[1] > y1) and (x[1] < y2)
#data structure model
class Information:
    """Contains general information"""
    def __init__(self): 
        self.elements_per_line = ELEMENTS_PER_LINE;
        self.number_of_elements = int(math.pow(self.elements_per_line,2))
        self.nodes_per_element = NODES_PER_ELEMENT
        self.number_of_nodes = int(math.pow(self.elements_per_line+1,2))
        self.rho = rho
        self.eps = eps

class Tables:
    """Contains all needed tables"""
    def __init__(self):
        #self.element_table = None
        self.local_to_global = None
        self.nodes_to_coordinates = None
        self.boundary_table = None

class Plots:
    """Contains plot data"""
    def __init__(self):
        self.xx = None
        self.yy = None
        self.zz = None
        self.x_max = ELEMENTS_PER_LINE * 2
        self.y_max = ELEMENTS_PER_LINE * 2
        self.step = STEP_SIZE  #visualisation mesh refinement

class Equations:
    """Contains equation data"""
    def __init__(self):
        self.K_h = None
        self.f_h = None
        self.u_h = None


class Data:
    """Container structure to be filled with data"""
    def __init__(self):
        self.info = Information()
        self.tables = Tables()
        self.equ = None
        self.plots = Plots()
        

def run():
    data = Data()
    
    with Timer('preparation'):
        #create element_table
        #data.tables.element_table = preparation.create_element_table(data.info) #unused?

        #import build
        import build_viertel as build
        #import build_ladung as build

        data.info.elements_per_line = -1


        data.info.number_of_elements = build.get_number_of_elements()
        data.info.nodes_per_element = build.get_nodes_per_element()
        data.info.number_of_nodes = build.get_number_of_nodes()

        #create local_to_global table
        #data.tables.local_to_global = preparation.create_local_to_global(data.info)
        data.tables.local_to_global = build.get_local_to_global().copy()
        
        #create nodes_to_coordinates
        #data.tables.nodes_to_coordinates = preparation.create_nodes_to_coordinates(data.info)
        data.tables.nodes_to_coordinates = build.get_nodes_to_coordinates().copy()
        print data.tables.nodes_to_coordinates


        #create boundary_table
        #data.tables.boundary_table = preparation.create_boundary_table(data.info)
        data.tables.boundary_table = build.get_boundary_table().copy()

        helper.write_matrix_to_file(data.tables.local_to_global, 'data/local_to_global.txt')
        helper.write_matrix_to_file(data.tables.nodes_to_coordinates, 'data/nodes_to_coordinates.txt')
        helper.write_matrix_to_file(data.tables.boundary_table, 'data/boundary_table.txt')


    with Timer('assembly'):
        #assemble equations
        data.equ = assembly.assemble_equations(data.info, data.tables)

        pass

    with Timer('solve'):
        #solve equations
        data.equ.u_h = solve.solve(data.equ)
        helper.write_matrix_to_file(data.equ.u_h, 'data/u_h.txt')
        print(data.equ.u_h)

    with Timer('plots'):
        #generate plot_data
        # data.plots.xx, data.plots.yy, data.plots.zz = plots.generate_plot_data(data.info, 
        #                             data.tables, data.plots, data.equ.u_h)

        #plots.contour_plot(data.info, data.tables, data.plots, data.equ.u_h)
        
        plots.contour_plot_new(data.info, 
                                     data.tables, data.plots, data.equ.u_h)

        #plots.fe_plot(data.info, data.tables, data.plots)
        pass


#execute run function when executed
if __name__ == '__main__':
    run()








