##Screw pinch convergence

*ConvergenceStudy.m* is a matlab script that will compute a finite difference estimate of the force gradient in a free boundary rotating ellipse for different values of perturbation delta R and compare it to SPEC force gradient. The maximal relative difference between both it plotted. To run, first be sure to have compiled **the debug version of SPEC** with

`>>> make dspec`

Leave the executable in SPEC source directory. Then, use the following command:

`>>> matlab -nodesktop < ConvergenceStudy.m`

Execution time is of the order of ~20 minutes. A plot *Convergence\_Torus.eps* is generated. You should observe a convergence till ~1E-10, finite difference errors are then forbidding a better agreement.





