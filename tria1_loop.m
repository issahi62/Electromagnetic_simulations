a =[]; 
B = zeros(100, 100); 
for t = 1:0.1:100 
for r = 1:100
B(:, r)= trial1(r, t);
end
end
figure(1);
imagesc(real(B))
%plot(real(a), 'LineWidth', 2); 

%% INSTANCE APPROACH

E = @(r,t)r.*exp(-1i*2.*t);

[R, T] = meshgrid(-1:0.01:1);
E = E(R, T); 
figure(2)
imagesc(real(E))
colormap
axis equal tight
title('E-field');
xlabel('x');
ylabel('y'); 