function write_rzdata(R,Z,outfname)

% Writes text file with R,Z data as two columns 
%
% INPUT
%   -R,Z      : data in the form of matrices (np x np)
%   -outfname : output file name with .txt form 
%
% written by J.Loizu (2016)

R1d=reshape(R,[1,size(R,1)*size(R,2)]);

Z1d=reshape(Z,[1,size(Z,1)*size(Z,2)]);

rz1d = [R1d; Z1d];

fileID = fopen(outfname,'w');

fprintf(fileID,'%1s %17s\n','R','Z');

fprintf(fileID,'%16.15f %16.15f\n',rz1d);

