function [ bns ] = fft_field( bn_raw, mpol)
% Transfer bn_raw to Fourier space
% mpol - number of poloidal harmonics

nfft = numel(bn_raw);
Y = fft(bn_raw,nfft);
P2 = -imag(Y/nfft);
P1 = P2(:,1:nfft/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);

bns = P1(1:mpol+1);


end

