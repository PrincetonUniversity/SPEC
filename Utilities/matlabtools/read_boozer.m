%% read_boozer( FILENAME, ROOT)
% =============================
%
% Reads hdf5 output from Boozer_xForms
% and produces data struct
%
% INPUT
% -----
%   -filename : name of Boz output 
%   -root     : root directory
%
% ------------------------------------%
% Written by S.Guinchard (05/22)      %
% from read_spec routine              %
% ------------------------------------%

function g_data = read_boozer(filename, root )
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
      g_data.([g_path{end}]) = read_boozer(filename, g_data_info.Groups(i).Name);
    end
  end
end
