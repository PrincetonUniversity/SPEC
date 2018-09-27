This constitutes a set of files that allow verification of free-boundary SPEC in a vacuum 5-field-period rotating ellipse field.

The file coils.fo.vacuum contains a segment representation (x,y,z) of a certain number of coils that generate the vacuum field. The current (fourth column) is interpreted in MA. 

The .sp file is the input for SPEC and contains the Rwc and Vns coefficients, which provide respectively the geometry of a control surface and the normal component of B on that surface. A Biot-Savart code was used to calculate B.n on the given control surface.

The solution obtained from SPEC shall be compared to the vacuum field solution obtained from Biot-Savart. One can compare the three components of B at different locations, mod(B) at different locations, Poincare plots, iota profiles, etc.

The two figures show a comparison between SPEC and Biot-Savart soltuions for the iota-profile and Poincare plot of the field lines. SPEc data plotted in black and Biot-Savart data in magenta or green. For the Poincare plot, the two data sets overlap almost perfectly so that to the eye the agreement is exact.

