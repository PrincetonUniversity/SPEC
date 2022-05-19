# Introduction

The class SPEC\_Namelist is a matlab tool to read, edit and write SPEC input files. Be careful though, always double check what kind of input file you generate - some errors are spotted by some tests in SPEC\_Namelist, but not everything is covered!

In this tutorial, we will show how to read a namelist, edit its Fourier resolution and change some specific Fourier harmonics, plot the initial guess and finally write the namelist.

# Reading a Namelist

## Adding SPEC\_Namelist path to MATLAB
First of all, make sure that the path to SPEC\_Namelist is provided to MATLAB. To do so, you can write
```Matlab
addpath(genpath('my/path/to/SPEC/Utilities/matlabtools'))
```

to check that MATLAB is using the right SPEC\_Namelist class, you can check which file it is using with the command
```Matlab
which SPEC_Namelist
```

## Read a SPEC input file
To create a SPEC\_Namelist instance from a SPEC input file, you need to give its absolute or relative path to the constructor of the class. For example, try 
```Matlab
nm = SPEC_Namelist('path/to/SPEC/InputFiles/TestCases/G3V08L3Fr.001.sp')
```

This will read the SPEC\_Namelist, create a MATLAB structure with all the input in a specific format, and check that some crucial inputs are provided. In case of missing inputs, it will print a few warnings:

```Matlab
nm = SPEC_Namelist('~/SPEC/InputFiles/TestCases/G3V08L3Fr.001.sp')
Warning: lboundary not provided. Setting with 0... 
> In SPEC_Namelist/initialize_structure (line 972)
  In SPEC_Namelist (line 73)
Warning: Missing Ivolume. Filling with 0... 
> In SPEC_Namelist/initialize_structure (line 1070)
  In SPEC_Namelist (line 73) 
Warning: Missing Isurf. Filling with 0... 
> In SPEC_Namelist/initialize_structure (line 1077)
  In SPEC_Namelist (line 73) 
Warning: Missing lboundary. Setting to zero 
> In SPEC_Namelist/initialize_structure (line 1122)
  In SPEC_Namelist (line 73) 
Warning: lboundary not provided. Setting with 0... 
> In SPEC_Namelist/initialize_structure (line 972)
  In SPEC_Namelist (line 73)
  In SPEC_Namelist/read_initial_guess (line 172)
  In SPEC_Namelist (line 81) 
Warning: Missing Ivolume. Filling with 0... 
> In SPEC_Namelist/initialize_structure (line 1070)
  In SPEC_Namelist (line 73)
  In SPEC_Namelist/read_initial_guess (line 172)
  In SPEC_Namelist (line 81) 
Warning: Missing Isurf. Filling with 0... 
> In SPEC_Namelist/initialize_structure (line 1077)
  In SPEC_Namelist (line 73)
  In SPEC_Namelist/read_initial_guess (line 172)
  In SPEC_Namelist (line 81) 
Warning: Missing lboundary. Setting to zero 
> In SPEC_Namelist/initialize_structure (line 1122)
  In SPEC_Namelist (line 73)
  In SPEC_Namelist/read_initial_guess (line 172)
  In SPEC_Namelist (line 81) 

nm = 

  SPEC_Namelist with properties:

              lists: {6x1 cell}
        physicslist: [1x1 struct]
        numericlist: [1x1 struct]
          locallist: [1x1 struct]
         globallist: [1x1 struct]
    diagnosticslist: [1x1 struct]
         screenlist: [1x1 struct]
      initial_guess: [1x1 struct]
```

## Data structure format
SPEC\_Namelist stores SPEC input data in nested structures. All structure and data names uses lower case - be careful whenever you access the data!

In addition, the class has an internal Fourier resolution, which is set by default to the largest poloidal and toroidal mode number in the input file. This resolution is a private attribute of the class and is not accessible to the user; to change it and increase or truncate all physical quantities in Fourier space to the resolution `mpol`, `ntor`, you can use

```Matlab
nm = nm.truncate_fourier_series( mpol, ntor )
```

All Fourier modes are stored in matrices of size `2*ntor+1 x mpol+1`. To get a specific mode, for example `Rbc(n,m)`, use
```Matlab
nm.get_fourier_harmonics( 'rbc', m, n )
```

# Modifying a SPEC\_Namelist
To modify any input, you can simply access it as you would with any MATLAB structure. For example, to change the toroidal flux in the third volume to one, write
```Matlab
nm.physicslist.tflux(3) = 1;
```

For the Fourier modes, setters routines are provided. For example, to change the mode Rbc(n,m) to a value `value`, write
```Matlab
nm = nm.set_fourier_harmonics( 'rbc', 1, 0, value );
```

If you want to set multiple modes, you can set them by providing arrays of values for `m`, `n` and `value`. Finally, if you want to use a computational boundary, a plasma boundary or the initial guess for the internal interfaces geometry from another SPEC input, you can use respectively `nm.set_boundary_from_namelist( filename, 'CB' )`, `nm.set_boundary_from_namelist( filename, 'PB' )` or `nm.read_initial_guess(filename)`. More information can be obtained by reading the documentation.

**DO NOT CHANGE HARMONICS WITHOUT THESE ROUTINES!** This might break the class and some of your input will be ignored.

# Visualizing the geometry of an input
Sometimes, SPEC crashes and it is difficult to know why. Some information can be obtained by plotting the geometry of computational boundary, the plasma boundary and the initial guess from the input file. 

## Plotting the plasma boundary
To plot the plasma boundary, write
```Matlab
nt = 1024; %Number of poloidal points
phi = 0;   %Toroidal angle
nm.plot_plasma_boundary( nt, phi, 1 );
```

The initial guess, if any, can be plotted with
```Matlab
nm.plot_initial_guess( nt, phi, 0 );
```

## Plotting the computational boundary
Similarly, the computational boundary can be plotted with
```Matlab
nm.plot_computational_boundary( nt, phi, 'N', 1 );
```

If you want to see the normal field on the computational boundary, switch `'N'` to (i) `'V'` to see the coils contribution, (ii) `'B'` to see the initial guess for the plasma contribution and (iii) `'F'` for the sum of both.


# Writing the namelist into a file
To write the Namelist, use the built in subroutine `write_input_file(filename)`. For example, do
```Matlab
nm.write_input_file( 'test.sp' )
```

will write an input file called test.sp with all the relevant quantities. Some input won\'t be written if not required - for example, if `nm.physicslist.Lfreebound` is set to zero, the fields `rwc`, `rws`, `zwc` and `zws` will not be written.









