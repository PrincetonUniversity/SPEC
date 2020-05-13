##Screw pinch convergence

*ConvergenceStudy.m* is a matlab script that will compute the force gradient in a free boundary plasma slab for different values of Lrad and compare it to the analytical force gradient. The maximal absolute difference between both it plotted. To run, first be sure to have compiled **the debug version of SPEC** with

`>>> make dspec`

Leave the executable in SPEC source directory. Then, use the following command:

`>>> matlab -nodesktop < ConvergenceStudy.m`

Execution time is of the order of ~5 minutes. A plot *Convergence\_Slab.eps* is generated. You should observe a convergence tille machine precision (~1E-15).





