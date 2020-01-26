function E = trial1(r, t)
%% E filed E(r,t)=E(r)exp(-iwot)
f = 1;
E_r = linspace(1, 2,100); 
w_0 = f;

E = E_r(r).*exp(1i*w_0.*t);

end