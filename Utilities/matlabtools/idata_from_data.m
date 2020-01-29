function idata = idata_from_data(data)

% Reads rotational transform data from field-line-tracing using output from SPEC
%
% INPUT
% - data  : data file obtained from read_spec(filename)
%
% OUTPUT
% - pdata : contains all the transform data, which can be fed into several routines for analyzing and ploting
%
% written by A. Baillod (2019)


nvol       = double(data.input.physics.Nvol);
Lfreebound = data.input.physics.Lfreebound;
nvol       = nvol + Lfreebound;

idata.iota = data.transform.fiota(:,2);
idata.sarr = data.transform.fiota(:,1);
idata.Mvol = nvol;
