function pdata = pdata_from_data(data)

% Reads Poincare data from field-line-tracing using output from SPEC
%
% INPUT
% - data   : data file obtained from read_spec(filename)
%
% OUTPUT
% - pdata  : contains all the poincare data, which can be fed into several routines for analyzing and ploting
%
% written by A.Baillod (2019)

nvol            = double(data.input.physics.Nvol);
Lfreebound      = data.input.physics.Lfreebound;
nvol            = nvol + Lfreebound;

pdata.R_lines   = data.poincare.R;
pdata.Z_lines   = data.poincare.Z;
pdata.npoinc    = data.input.diagnostics.nPpts;
pdata.th_lines  = data.poincare.t;
pdata.rho_lines = data.poincare.rho;
pdata.Igeometry = data.input.physics.Igeometry;
pdata.Lfreebound= data.input.physics.Lfreebound;
pdata.mn        = data.output.mn;
pdata.in        = data.output.in;
pdata.im        = data.output.im;
pdata.Rbc       = data.output.Rbc;
pdata.Rbs       = data.output.Rbs;
pdata.Zbc       = data.output.Zbc;
pdata.Zbs       = data.output.Zbs;
pdata.Nfp       = data.input.physics.Nfp;
pdata.Mvol      = nvol;
pdata.Nvol      = nvol - Lfreebound;

