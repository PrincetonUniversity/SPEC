% Interpolate N-D array with arbitrary mixture of upsampling and downsampling
%   using Fourier interpolation. Downsampling symmetrically truncates outer 
%   portions of k-space. Upsampling uses zero-filling.
%
%   Usage : out = fftInterpolate(in,newsz)
%
%       out   : output N-dimensional data with size(out) == newsz
%
%       in    : input N-dimensional data
%       newsz : desired interpolated size of output data
%
%   Example :
%
%       fftInterpolate(phantom(32),[24 40]) interpolates a 32x32 image to 24x40         
% 
%   Copyright 2008 Matthias Christian Schabel (matthias @ stanfordalumni . org)
%   University of Utah Department of Radiology
%   Utah Center for Advanced Imaging Research
%   729 Arapeen Drive
%   Salt Lake City, UT 84108-1218
%   
%   2010/12/10 MCS - fixed and simplified

function out = fftInterpolate(in,newsz)

sz = size(in);
nd = ndims(in);

% no interpolation
if (newsz == sz) out = in; return; end;

% negative or zero scaling factor is meaningless
if (any(newsz <= 0)) error('fftInterpolate :: bad size'); end;

% do it in one fell swoop
fac = prod(newsz./sz);

f = find(newsz > sz);

center = floor(sz/2)+1;
lo = floor(center-newsz/2);
hi = lo+newsz-1;

lo(f) = 1;
hi(f) = sz(f);

rng = [lo; hi]';

inft = subRange(fftshift(fftn(in)),rng);

out = zeros(newsz);

sz = size(inft);

% downscaling
centerd = floor(newsz/2);
lod = ceil(centerd-newsz/2)+1;
hid = lod+newsz-1;

% upscaling
centeru = floor(newsz/2)+1+odd(newsz);
lou = ceil(centeru-sz/2);
hiu = lou+sz-1;

lo = lod;
hi = hid;

% merge down and upscaling to get final range
lo(f) = lou(f);
hi(f) = hiu(f);

rng = [lo; hi]';

out = fac*ifftn(fftshift(assignSubRange(out,rng,inft)));

% eliminate residual complex values if input array is real
if (isreal(in))
    out = real(out);
end;

return;

function outodd = odd(inodd)
  outodd = mod(inodd,2);
return;

function subData = subRange(data,rng)

if (isempty(rng)) 
    subData = data; 
    return; 
end;

sz = size(data);
dim = length(sz);

lo = rng(:,1)';
hi = rng(:,2)';

if (length(lo) ~= dim || length(hi) ~= dim)
    error('subRange :: dimension mismatch');
end;

% replace zeros with lower/upper limit 
for i=1:dim
    if (lo(i) == 0) lo(i) = 1; end;
    if (hi(i) == 0) hi(i) = sz(i); end;
end;

if (any(lo<1) | any(hi > sz))
    error('subRange :: sub-range out of bounds');
end;

switch (dim)
    case 1,  subData = data(lo(1):hi(1)); 
    case 2,  subData = data(lo(1):hi(1),...
                           lo(2):hi(2)); 
    case 3,  subData = data(lo(1):hi(1),...
                           lo(2):hi(2),...
                           lo(3):hi(3)); 
    case 4,  subData = data(lo(1):hi(1),...
                           lo(2):hi(2),...
                           lo(3):hi(3),...
                           lo(4):hi(4)); 
    case 5,  subData = data(lo(1):hi(1),...
                           lo(2):hi(2),...[]
                           lo(3):hi(3),...
                           lo(4):hi(4),...
                           lo(5):hi(5)); 
    otherwise,
            % generate string and use eval
            str = 'subData = data(';
            for i=1:dim
                str = [str 'lo(' num2str(i) '):hi(' num2str(i) ')'];
                if (i~=dim)
                    str = [str ','];
                else
                    str = [str ');'];
                end;
            end;
            eval(str);
end;
        
return;


function out = assignSubRange(out,rng,in)

if (isempty(rng)) return; end;

sz = size(out);
dim = length(sz);

lo = rng(:,1)';
hi = rng(:,2)';

if (length(lo) ~= dim || length(hi) ~= dim)
    error('assignSubRange :: dimension mismatch');
end;

% replace zeros with lower/upper limit 
for i=1:dim
    if (lo(i) == 0) lo(i) = 1; end;
    if (hi(i) == 0) hi(i) = sz(i); end;
end;

if (any(lo<1) | any(hi > sz))
    error('assignSubRange :: sub-range out of bounds');
end;

if (size(in) ~= (hi-lo+1))
    error('assignSubRange :: input data size mismatch');
end;

switch (dim)
    case 1,  out(lo(1):hi(1)) = in; 
    case 2,  out(lo(1):hi(1),...
                 lo(2):hi(2)) = in; 
    case 3,  out(lo(1):hi(1),...
                 lo(2):hi(2),...
                 lo(3):hi(3)) = in; 
    case 4,  out(lo(1):hi(1),...
                 lo(2):hi(2),...
                 lo(3):hi(3),...
                 lo(4):hi(4)) = in; 
    case 5,  out(lo(1):hi(1),...
                 lo(2):hi(2),...
                 lo(3):hi(3),...
                 lo(4):hi(4),...
                 lo(5):hi(5)) = in; 
    otherwise,
            % generate string and use eval
            str = 'out(';
            for i=1:dim
                str = [str 'lo(' num2str(i) '):hi(' num2str(i) ')'];
                if (i~=dim)
                    str = [str ','];
                else
                    str = [str ') = in;'];
                end;
            end;
            eval(str);
end;
        
return;


