import numpy as np

def solve(equ):
    """
    Solves the equation system K_h * u_h = f_h.
    """
    return np.linalg.solve(equ.K_h, equ.f_h)



#execute run function when executed
if __name__ == '__main__':
    import fem
    fem.run()
