function d = Dfac(tf,pf,mu,del)

mub   = mu*del/2;
sigma = sqrt(abs(mu^2-1));
sigmb = sigma*del/2;

if(mu>=1)

 d = (mu^2/(32*pi^2))*sigma*(pf^2*cot(mub)^2-tf^2)*(tan(sigmb)+cot(sigmb));
    
elseif(mu<1)
 
 if(mu==0)
  d = (pf^2/(4*pi^2*del^2))*(sigma/sinh(2*sigmb));
 else
  d = (mu^2/(32*pi^2))*(2*sigma/sinh(2*sigmb))*(pf^2*cot(mub)^2-tf^2);
 end

end

