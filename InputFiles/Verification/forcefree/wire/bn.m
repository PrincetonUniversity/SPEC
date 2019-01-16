function [ theta, bn_raw ] = bn( Rbc, Zbs, n, R0, mu0I )
% calculate the Fouier transform of Bn at given boundary
% Rbc  - R harmonics of theta
% Zbs  - Z harmonics of theta
% n    - number of points
% R0   - major radius of the wire
% mu0I - mu0*current

nfft = n;

theta = linspace(0, 2*pi, nfft);

Rij = zeros(size(theta));
Zij = Rij;
dRij = Rij;
dZij = Zij;

bn_raw = zeros(size(linspace(0, 2*pi, nfft)));

for i = 1: numel(Rbc)
    Rij = Rij + Rbc(i) * cos((i-1) * theta);
    Zij = Zij + Zbs(i) * sin((i-1) * theta);
    dRij = dRij - Rbc(i) * (i-1) * sin((i-1) * theta);
    dZij = dZij + Zbs(i) * (i-1) * cos((i-1) * theta);
end

for i = 1: nfft
    [Br, Bz] = field( R0, mu0I, Rij(i), Zij(i) );
    %dr = Rij(i+1) - Rij(i);
    %dz = Zij(i+1) - Zij(i);
    dr = dRij(i);
    dz = dZij(i);
    %dd = sqrt(dr^2 + dz^2);
    nr = dz*Rij(i); %/ dd;
    nz = -dr*Rij(i); %/ dd;
    bn_raw(i) = Br*nr + Bz * nz;
end


%Y = fft(bn_raw,nfft);
%P2 = abs(Y/nfft);
%P1 = P2(:,1:nfft/2+1);
%P1(:,2:end-1) = 2*P1(:,2:end-1);

%bns = P1(1:mpol);

end

