function fdata = fdata_from_data(data)

% Reads Poincare data from field-line-tracing using output from SPEC
%
% INPUT
% - data : data file obtained from read_spec(filename)
%
% OUTPUT
% - pdata    : contains all the poincare data, which can be fed into several routines for analyzing and ploting
%
% written by J.Loizu (2017)

global machform;
machform = 's'; 

nvol     = double(data.input.physics.Nvol);
Lfreebound = data.input.physics.Lfreebound;

fdata    = data.input.physics;

f = fieldnames(data.output);
for i = 1:length(f)
fdata.(f{i}) = data.output.(f{i});
end

f = fieldnames(data.vector_potential);
for i = 1:length(f)
fdata.(f{i}) = data.vector_potential.(f{i});
end

fdata.Mregular = data.input.numerics.Mregular;