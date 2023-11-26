
addpath('airfoils');
FOIL=importselig('n2415.dat');
foil = xfoilCl(FOIL.normalized,0.346,1e6,0)
