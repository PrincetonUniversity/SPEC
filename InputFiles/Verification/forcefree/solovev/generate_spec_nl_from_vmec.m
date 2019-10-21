function generate_spec_nl_from_vmec(output_name, vmec_file, spec_template, Nvol, tflux, imode)
% Build a SPEC namelist based on a VMEC output wout file
% Currenlt, only for stellarator symmetry case,
% only zero current is supported. This script will only generate
% an initial guess. One needs to use run_spec_iter_fixcurr to find solution
% with zero net current.
%
% Prerequisite: matlabVMEC (v2.92) by Samuel Lazerson
% https://www.mathworks.com/matlabcentral/fileexchange/29031-matlabvmec
%
% Inputs:
% output_name - the name of the output namelist
% vmec_file - the path of VMEC wout file
% spec_template - the path of a SPEC template
% Nvol - number of volumes
% tflux - an array of tflux enclosed by each interface
% imode - 0:zero current iteration
%         1:fix rotational transform
%
% Zhisong Qu 3/12/2018

% default for namelist
Lrad = 8;
nPtrj = -1;

% some numerical settings for integration
nint = 1001;

f = read_vmec(vmec_file);

ntor = f.ntor
mpol = f.mpol
mn = f.mnmax
im = f.xm;
in = f.xn/f.nfp;

%dtflux = diff(tflux);
phiedge = f.phi(end);
normphi = f.phi / phiedge;

% obtain pressure and iota by intepolation
pressure = zeros([1,Nvol]);
iotasum = zeros([1,Nvol]);
for i = 1:Nvol
    if i==1
        phistart = 0.0;
    else
        phistart = tflux(i-1);
    end
    phiend = tflux(i);
    
    phiint = linspace(phistart, phiend, nint);
    pressureint = trapz(phiint, interp1(normphi, f.presf, phiint, 'pchib'));
    iotaint = trapz(phiint, interp1(normphi, f.iotaf, phiint, 'pchib'));
    pressure(i) = pressureint ./ (phiend - phistart);
    iotasum(i) = iotaint;
end

% multiply by mu
pressure = pressure * 4 * pi * 1e-7;
% normalize
p0 = pressure(1);
if (p0 < 1e-10)
    pressure = pressure * 0;
else
    pressure = pressure / pressure(1);
end

% calculate pflux
pflux = zeros(size(tflux));
pflux(1) = 0;
for i = 2:Nvol
    pflux(i) = -sum(iotasum(2:i));
end

% construct rac, zas
rac = f.rmnc(1:ntor+1,1);
zas = f.zmns(1:ntor+1,1);

% construct rbc, zbc
% we need to flip poloidal angle
rbc = zeros([1,mn]);
zbs = zeros([1,mn]);
for i = 1:mn
    if im(i) == 0
        rbc(i) = f.rmnc(i,end);
        zbs(i) = f.zmns(i,end);
    else
        rbc(i) = f.rmnc(i,end);
        zbs(i) = -f.zmns(i,end);
    end
end

% construct inner surfaces
sizerbc = size(f.rmnc);
sizerbc(2) = Nvol+1;
irbc = zeros(sizerbc);
izbs = zeros(sizerbc);

irbc(:,1) = f.rmnc(:,1);
izbs(:,1) = f.zmns(:,1);
irbc(:,end) = rbc;
izbs(:,end) = zbs;

for i = 1:sizerbc(1)
    if im(i) == 0
        irbc(i,2:end-1) = interp1(normphi, f.rmnc(i,:), tflux(1:end-1));
        izbs(i,2:end-1) = interp1(normphi, f.zmns(i,:), tflux(1:end-1));
    else
        irbc(i,2:end-1) = interp1(normphi, f.rmnc(i,:), tflux(1:end-1));
        izbs(i,2:end-1) = -interp1(normphi, f.zmns(i,:), tflux(1:end-1));
    end
end

% finially, we need to filp the sign of n for m = 0
in(1:ntor+1) = abs(in(1:ntor+1));

specnl.sref{1}.name = 'pressure';
specnl.sref{1}.data = pressure;
specnl.sref{2}.name = 'tflux';
specnl.sref{2}.data = tflux;
specnl.sref{3}.name = 'pflux';
specnl.sref{3}.data = pflux;
specnl.sref{4}.name = 'Lfindzero';
specnl.sref{4}.data = 2;
specnl.sref{5}.name = 'Lrad';
specnl.sref{5}.data = Lrad * ones(size(pflux), 'int32');
specnl.sref{6}.name = 'nPtrj';
specnl.sref{6}.data = -1 * ones(size(pflux), 'int32');
specnl.sref{7}.name = 'Nvol';
specnl.sref{7}.data = Nvol;
specnl.sref{8}.name = 'phiedge';
specnl.sref{8}.data = phiedge;
specnl.sref{9}.name = 'pscale';
specnl.sref{9}.data = p0;
specnl.sref{10}.name = 'Mpol';
specnl.sref{10}.data = mpol;
specnl.sref{11}.name = 'Ntor';
specnl.sref{11}.data = ntor;
specnl.sref{12}.name = 'Linitialize';
specnl.sref{12}.data = 0;

specnl.sref{13}.name = 'Rac';
specnl.sref{13}.data = rac;
specnl.sref{14}.name = 'Zas';
specnl.sref{14}.data = zas;
specnl.sref{15}.name = 'Ras';
specnl.sref{15}.data = zeros(size(rac), 'double');
specnl.sref{16}.name = 'Zac';
specnl.sref{16}.data = zeros(size(rac), 'double');

specnl.sref{17}.name = 'Lconstraint';
specnl.sref{17}.data = imode;

% specify mu
specnl.sref{18}.name = 'mu';
if imode == 0
    specnl.sref{18}.data = zeros(size(pflux), 'double');
else
    specnl.sref{18}.data = 0.1 * ones(size(pflux), 'double');
end
% specify rotational transform
if imode == 1
    specnl.sref{19}.name = 'iota';
    specnl.sref{20}.name = 'oita';
    
    iota = zeros([1,Nvol+1]);
    iota(2:end) = interp1(normphi, f.iotaf, tflux(1:end));
    
    specnl.sref{19}.data = -iota;
    specnl.sref{20}.data = -iota;
end

specnl.sref2d{1}.name = 'Rbc';
specnl.sref2d{1}.data = rbc;
specnl.sref2d{2}.name = 'Zbs';
specnl.sref2d{2}.data = zbs;
specnl.sref2d{3}.name = 'Rbs';
specnl.sref2d{3}.data = zeros(size(rbc), 'double');
specnl.sref2d{4}.name = 'Zbc';
specnl.sref2d{4}.data = zeros(size(rbc), 'double');

specnl.interfaces.irbc = irbc;
specnl.interfaces.izbs = izbs;
specnl.interfaces.irbs = zeros(size(irbc), 'double');
specnl.interfaces.izbc = zeros(size(irbc), 'double');
specnl.interfaces.im = im;
specnl.interfaces.in = in;
specnl.interfaces.mn = mn;
specnl.interfaces.nvol = Nvol;

write_spec_input_vmec(spec_template, output_name, specnl);
end

