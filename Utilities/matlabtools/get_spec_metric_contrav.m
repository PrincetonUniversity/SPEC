function gmatcon = get_spec_metric_contrav(data,lvol,sarr,tarr,zarr)
 
% GET_SPEC_METRIC_CONTRAV( DATA, LVOL, SARR, TARR, ZARR )
% =======================================================
%
% Calculates contravariant metric elements in volume number lvol
%
% INPUT
% -----
%   -data   must be produced by calling e.g. read_spec(filename)
%   -lvol   volume number
%   -sarr   is the array of values for the s-coordinate
%   -tarr   is the array of values for the theta-coordinate
%   -zarr   is the array of values for the zeta-coordinate
%
% OUTPUT
% ------
%   gmatcon{i,j}(s,theta,zeta)
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2017)

Istellsym = data.input.physics.Istellsym;
if Istellsym~=1
    error('Only valid for stellarator symmetric equilibria')
end


% Auxiliary variables

ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);

% First get the covariant metric matrix
gmatcov = get_spec_metric(data,lvol,sarr,tarr,zarr);

% Allocate space for contravariant metric matrix
gmatcon = cell(3,3);  

% Initiate calculation loop (over each position on the surface)

for is=1:ns
 for it=1:nt
  for iz=1:nz


  % Define arrays of gcov, gcon, and rhs

  rhs     = transpose([1 1 1 0 0 0]);


  % Define transformation matrix G
  G       = zeros(6,6);

  G(1,1)  = gmatcov{1,1}(is,it,iz);
  G(2,2)  = gmatcov{2,2}(is,it,iz);
  G(3,3)  = gmatcov{3,3}(is,it,iz);
  G(4,4)  = gmatcov{1,1}(is,it,iz);
  G(5,5)  = gmatcov{1,1}(is,it,iz);
  G(6,6)  = gmatcov{2,2}(is,it,iz);

  G(1,4)  = gmatcov{1,2}(is,it,iz);
  G(1,5)  = gmatcov{1,3}(is,it,iz);

  G(2,4)  = gmatcov{1,2}(is,it,iz);
  G(2,6)  = gmatcov{2,3}(is,it,iz);
 
  G(3,5)  = gmatcov{1,3}(is,it,iz);
  G(3,6)  = gmatcov{2,3}(is,it,iz);

  G(4,2)  = gmatcov{1,2}(is,it,iz);
  G(4,6)  = gmatcov{1,3}(is,it,iz);

  G(5,3)  = gmatcov{1,3}(is,it,iz);
  G(5,6)  = gmatcov{1,2}(is,it,iz);

  G(6,3)  = gmatcov{2,3}(is,it,iz);
  G(6,5)  = gmatcov{1,2}(is,it,iz);

  % Calculate gcon elements
  gcon    = G\rhs; 

  % Construct contravariant metric matrix

  gmatcon{1,1}(is,it,iz) = gcon(1);  %g^ss
  gmatcon{2,2}(is,it,iz) = gcon(2);  %g^tt
  gmatcon{3,3}(is,it,iz) = gcon(3);  %g^zz
  gmatcon{1,2}(is,it,iz) = gcon(4);  %g^st
  gmatcon{1,3}(is,it,iz) = gcon(5);  %g^sz
  gmatcon{2,3}(is,it,iz) = gcon(6);  %g^tz
  
  gmatcon{2,1}(is,it,iz) = gmatcon{1,2}(is,it,iz);  % by symmetry of g
  gmatcon{3,1}(is,it,iz) = gmatcon{1,3}(is,it,iz);
  gmatcon{3,2}(is,it,iz) = gmatcon{2,3}(is,it,iz);

  end
 end
end

