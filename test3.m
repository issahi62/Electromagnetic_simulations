function A = test3(x,y)
r = sqrt(x.^2 + y.^2);
k = 2 * pi;
w_0 = 1;
t = linspace(0,1,length(x));
A = exp(1i*k.*r - w_0.*t);
end