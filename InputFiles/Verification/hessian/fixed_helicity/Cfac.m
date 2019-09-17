function c = Cfac(pm,tf,pf,mu,del)

mub   = mu*del/2;
sigma = sqrt(abs(mu^2-1));
sigmb = sigma*del/2;

if(mu>=1)

 c = (mu^2/(32*pi^2))*( pm*sigma*(tan(sigmb)-cot(sigmb))*(tf + pm*pf*cot(mub))^2 + 2*mu*(pm*(tf^2-pf^2)*cot(mub) + tf*pf*(cot(mub)^2-1)) );
    
elseif(mu<1)

 if(mu==0)
  c = -pm*sigma*coth(2*sigmb)*pf^2/(4*pi^2*del^2);
 else
  c = (mu^2/(32*pi^2))*( -pm*2*sigma*coth(2*sigmb)*(tf + pm*pf*cot(mub))^2 + 2*mu*(pm*(tf^2-pf^2)*cot(mub) + tf*pf*(cot(mub)^2-1)) );
 end
end
