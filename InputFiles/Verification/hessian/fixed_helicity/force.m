function f = force(x,R3,mu,tf,pf,bcont)

mub(1)  = mu(1)*x(1)/2;
mub(2)  = mu(2)*x(2)/2;
mub(3)  = mu(3)*(R3-x(1)-x(2))/2;

if(bcont==1)
 pf(1) = -mu(2)*tf(2)*0.5*x(1);
 pf(3) = mu(2)*tf(2)*0.5*(R3-x(1)-x(2));
end

gmu  = zeros(1,3);
del  = [x(1) x(2) R3-x(1)-x(2)];

for i=1:3
 if(mu(i)==0)
  gmu(i) = (2/del(i))^2;
 else
  gmu(i) = mu(i)^2 / sin(mub(i))^2;
 end
end

f(1) = gmu(2) * (tf(2)^2 + pf(2)^2) - gmu(1) * (tf(1)^2 + pf(1)^2) ;

f(2) = gmu(3) * (tf(3)^2 + pf(3)^2) - gmu(2) * (tf(2)^2 + pf(2)^2) ;

%f(1) = mu(2)^2 * (tf(2)^2 + pf(2)^2) / sin(mub(2))^2 - mu(1)^2 * (tf(1)^2 + pf(1)^2) / sin(mub(1))^2 ;

%f(2) = mu(3)^2 * (tf(3)^2 + pf(3)^2) / sin(mub(3))^2 - mu(2)^2 * (tf(2)^2 + pf(2)^2) / sin(mub(2))^2 ;
