function data = read_spec(filename)
  
% Reads the HDF5 output file produced by SPEC
%
% INPUT
% -  filename : path to the HDF5 output file (e.g. 'testcase.h5')
%
% OUTPUT
% -  data     : contains all data from the SPEC run, which can be fed into several routines for analyzing and plotting
%
% written by J.Schilling (2019)
% based on the routines by J.Loizu and read_hdf5.m by S. Lazerson

% Try to read the file first
try
    data_info = h5info(filename,'/');
catch h5info_error
    data=-1;
    disp(['ERROR: Opening HDF5 File: ' filename]);
    disp(['  -identifier: ' h5info_error.identifier]);
    disp(['  -message:    ' h5info_error.message]);
    disp('      For information type:  help read_hdf5');
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

end
