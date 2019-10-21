function h = helicity(mub,tf,pf)

if(mub==0)
 h = 0;
else
 h = (tf^2 + pf^2)*mub/(2*sin(mub)^2);
end
