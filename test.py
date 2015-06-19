#      
# Copyright (c) 2015, Peishi Jiang

import GreenAmpt

def main():
    """
    test
    """
    K = 1.09; psi = 11.01; dtheta = 0.247; dt = 0.166
    i = [1.08, 1.26, 1.56, 1.92, 2.22, 2.58, 3.84, 6.84, 19.08, \
         9.90, 4.86, 3.12, 2.52, 2.16, 1.68, 1.44, 1.14, 1.02]
    a = GreenAmpt(K, dt, dtheta, psi, i)
    F, f, t =a.F_f()    

if __name__ == "__main__":
    main()