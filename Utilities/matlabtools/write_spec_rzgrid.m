function write_spec_rzgrid(data, nz0,lvol,outfname)


%
% WRITE_SPEC_RZGRID( DATA, NZ0, LVOL, OUTFNAME )
% ==============================================
%
% Writes text file with grid coordinate data points R,Z on a given volume, 
% on a toroidal plane, as two columns 
%
% INPUT
% -----
%   -fname    : input file name with .h5 form
%   -nz0      : toroidal plane number
%   -lvol     : volume number
%   -outfname : output file name with .txt form 
%
% written by J.Loizu (2018)


Nt    = data.grid.Nt;
Nz    = data.grid.Nz;
iz    = nz0-1;
Rij   = data.grid.Rij;
Zij   = data.grid.Zij;

Rdata = Rij(lvol,1+Nt*iz:(iz+1)*Nt,:);
Zdata = Zij(lvol,1+Nt*iz:(iz+1)*Nt,:);

R1d   = reshape(Rdata,[1,size(Rdata,2)*size(Rdata,3)]);

Z1d   = reshape(Zdata,[1,size(Zdata,2)*size(Zdata,3)]);

rz1d = [R1d; Z1d];

fileID = fopen(outfname,'w');

fprintf(fileID,'%1s %17s\n','R','Z');

fprintf(fileID,'%16.15f %16.15f\n',rz1d);

