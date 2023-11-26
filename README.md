This folder accumulates the required functions and XFOIL.exe file required to run the optimization
routine of a wing. This folder should be added to the MATLAB path before Wing_optimize.m is exectuted

Wing_optimize is the header file which all will execute the entire optimization routine.

the three design variables are shown below and are formed into the array which is optimized iteratively:


design variables
   b    Cr    Ct 
_______________________
[ x(1), x(2), x(3) ]


The output of the file is within the command line.

the first three values printed in the large array represent the design variables for wing geometry
the last three represent optimization output parameters and convergence detection.


The figures include the geometric shape of the selected airfoil
the values of each constraint as it iterated through the optimization

Figure 3 displays the value of the magnitude of the force of drag (minimization parameter) as it iterates 
through the optimization routine. 

The starting values and constraints that are saved here result in convergence on a solution after 
289 iterations.

