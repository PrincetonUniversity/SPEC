function data = read_spec(filename)

% 
% READ_SPEC ( FILENAME )
% ======================
%
% Reads the HDF5 output file produced by SPEC
%
% INPUT
% -----
% -  filename : path to the HDF5 output file (e.g. 'testcase.h5')
%
% OUTPUT
% ------
% -  data     : contains all data from the SPEC run, which can be fed into several routines for analyzing and plotting
%
% written by J.Schilling (2019)
% based on the read_spec_<quantity>.m routines by J. Loizu and read_hdf5.m by S. Lazerson

% Try to read the file first
try
    h5info(filename,'/');
catch h5info_error
    data=-1;
    disp(['ERROR: Opening HDF5 File: ' filename]);
    disp(['  -identifier: ' h5info_error.identifier]);
    disp(['  -message:    ' h5info_error.message]);
    return
end

% recursive reading routine
function g_data = getGroup(filename, root)
  g_data_info = h5info(filename, root);
  ngroups     = length(g_data_info.Groups);
  nvars       = length(g_data_info.Datasets);
  % Get datasets in root node
  for i = 1: nvars
    g_data.([g_data_info.Datasets(i).Name]) = h5read(filename,[root '/' g_data_info.Datasets(i).Name]);
    natts = length(g_data_info.Datasets(i).Attributes);
    for j=1:natts
        g_data.([g_data_info.Datasets(i).Attributes(j).Name]) = g_data_info.Datasets(i).Attributes(j).Value{1};
    end
  end
  % get groups in root node
  if ngroups > 0
    for i = 1 : ngroups
      g_path = strsplit(g_data_info.Groups(i).Name, '/');
      g_data.([g_path{end}]) = getGroup(filename, g_data_info.Groups(i).Name);
    end
  end
end

% start recursion at root node
data = getGroup(filename, '/');

% make adjustments for compatibility with previous reading routines
Nvol = data.input.physics.Nvol;
Mvol = data.output.Mvol;
Lrad = data.input.physics.Lrad;

% vector potential
cAte = cell(Nvol, 1);
cAto = cell(Nvol, 1);
cAze = cell(Nvol, 1);
cAzo = cell(Nvol, 1);

% grid
cRij = cell(Nvol, 1);
cZij = cell(Nvol, 1);
csg  = cell(Nvol, 1);
cBR  = cell(Nvol, 1);
cBp  = cell(Nvol, 1);
cBZ  = cell(Nvol, 1);

% Poincare data
cT   = cell(Nvol, 1);
cRho = cell(Nvol, 1);
cR   = cell(Nvol, 1);
cZ   = cell(Nvol, 1);

% split into separate cells for nested volumes
start=1;
for i=1:Mvol
  % vector potential
  cAte{i} = data.vector_potential.Ate(start:start+Lrad(i),:);
  cAto{i} = data.vector_potential.Ato(start:start+Lrad(i),:);
  cAze{i} = data.vector_potential.Aze(start:start+Lrad(i),:);
  cAzo{i} = data.vector_potential.Azo(start:start+Lrad(i),:);

  % grid
  cRij{i} = data.grid.Rij(start:start+Lrad(i),:)';
  cZij{i} = data.grid.Zij(start:start+Lrad(i),:)';
  csg{i}  = data.grid.sg(start:start+Lrad(i),:)';
  cBR{i}  = data.grid.BR(start:start+Lrad(i),:)';
  cBp{i}  = data.grid.Bp(start:start+Lrad(i),:)';
  cBZ{i}  = data.grid.BZ(start:start+Lrad(i),:)';
 
  % move along combined array dimension
  start = start + Lrad(i)+1;
end

% check if any poincare data was written to the output file
if isfield(data, 'poincare')
  % copy to get dimensions right
  data.poincare.rho = data.poincare.s;

  % disentangle Poincare data and generate rho entry; all this code is necessary since it depends on the volume index...
  start=1;
  for i=1:Nvol
    nPtrj = data.input.diagnostics.nPtrj(i);
    if (nPtrj==-1)
      nPtrj = 2*Lrad(i);
    end

    % In the innermost volume, there are exactly as many trajectories as specified.
    if ~(data.input.physics.Igeometry==1 || i>1)
      % mimic LREGION() macro functionality and the corresponding logic in pp00aa
      nPtrj = nPtrj-1;
    end

    % disp(sprintf('%d: %d ... %d', i, start, start+nPtrj));

    % create rho entry for Poincare plot
    offset = double(i-1)./double(Nvol);
    data.poincare.rho(:,:,start:start+nPtrj) = 0.5*(data.poincare.s(:,:,start:start+nPtrj)+1.0)./double(Nvol)+offset;
  
    % In all the outer volumes (for Nvol>1), there is one additional trajectory than specified in the input file.
    start = start + nPtrj+1;
  end

  % ensure compatibility with Joaquim's former reading routine
  data.poincare.t   = permute(data.poincare.t,   [3,1,2]);
  data.poincare.rho = permute(data.poincare.rho, [3,1,2]);
  data.poincare.R   = permute(data.poincare.R,   [3,1,2]);
  data.poincare.Z   = permute(data.poincare.Z,   [3,1,2]);
end % check presence of data.poincare


% replace original content in data structure
data.vector_potential.Ate = cAte;
data.vector_potential.Ato = cAto;
data.vector_potential.Aze = cAze;
data.vector_potential.Azo = cAzo;

data.grid.Rij = cRij;
data.grid.Zij = cZij;
data.grid.sg  = csg;
data.grid.BR  = cBR;
data.grid.Bp  = cBp;
data.grid.BZ  = cBZ;


end

