function [outdata] = resample2D(indata,Ntnew,Nznew)
%Resample the data using FFT

outdata = fftInterpolate(indata, [Ntnew, Nznew]);

%[Ntold, Nzold] = size(indata);
%fftold = fft2(indata);
%fftout = fftshift(fftold);

%fftnew = zeros(Ntnew, Nznew);
%lt = min(Ntnew/2, Ntold/2);
%lz = min(Nznew/2, Nzold/2);

%fftnew(Ntnew/2-lt+1:Ntnew/2+lt, Nznew/2-lz+1:Nznew/2+lz) = fftout(Ntold/2-lt+1:Ntold/2+lt, Nzold/2 - lz+1:Nzold/2 + lz);
%fftnew = ifftshift(fftnew);

%fftnew = zeros(Ntnew, Nznew);

%fftnew(1:min(Ntnew/2, Ntold/2), 1:min(Nznew/2, Nzold/2)) = fftout(1:min(Ntnew/2, Ntold/2), 1:min(Nznew/2, Nzold/2));
%fftnew(1:min(Ntnew/2, Ntold/2), Nznew - min(Nznew/2, Nzold/2)+1:end) = fftout(1:min(Ntnew/2, Ntold/2), Nzold - min(Nznew/2, Nzold/2)+1:end);
%fftnew(Ntnew - min(Ntnew/2, Ntold/2)+1:end, 1:min(Nznew/2, Nzold/2)) = fftout(Ntold - min(Ntnew/2, Ntold/2)+1:end, 1:min(Nznew/2, Nzold/2));
%fftnew(Ntnew - min(Ntnew/2, Ntold/2)+1:end, Nznew - min(Nznew/2, Nzold/2)+1:end) = fftout(Ntold - min(Ntnew/2, Ntold/2)+1:end, Nzold - min(Nznew/2, Nzold/2)+1:end);

%fftnew = fftnew / Ntold / Nzold * Ntnew * Nznew;

%outdata = real(ifft2(fftnew));
end

