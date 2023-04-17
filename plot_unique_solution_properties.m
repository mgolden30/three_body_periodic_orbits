%{
Make plots of unique solutions
%}

unique = [1,2,4,5,7,8,9,22,23,33,67,90,91,106,108,156];


for i = 1:numel(unique)

load( "solutions/" + unique(i) + ".mat");

figure(1);
plot_state(y,1);
set( gcf, 'Color' ,'w');
drawnow
saveas(gcf, "unique_solutions/traj_" + unique(i) + ".png");



figure(3);
scatter( 1:N, abs(1 + hamiltonian( reshape(y(1:8*N), [8,N])) ), 'filled' );
set(gca, 'yscale', 'log');
xlabel('timestep');
ylabel('|H+1|');
set( gcf, 'Color' ,'w');
saveas(gcf, "unique_solutions/delta_H_" + unique(i) + ".png");



% Fourier spectra of y
figure(4);
spectra = fft( reshape(y(1:8*N), [8,N] ), [], 2 )/N;
nice = abs(spectra);
xx = 1:(N/2-1);
nice = nice(:,xx);
scatter( xx-1, max(nice), 'filled' );
set( gca, 'YScale', 'log')
ylabel('max|{\bf z}_\omega|');
xlabel('\omega')
xlim([0, N/2-1]);
drawnow
set( gcf, 'Color' ,'w');
saveas(gcf, "unique_solutions/fourier_coeffs_" + unique(i) + ".png");


% angular momentum 
x = reshape( y(1:8*N), [8,N] ); 
r1 = x(1:2, :);
r2 = x(3:4, :);
p1 = x(5:6, :);
p2 = x(7:8, :);
r3 = -r1-r2;
p3 = -p1-p2;
cross = @(x,y) x(1,:).*y(2,:) - x(2,:).*y(1,:);

L = cross(r1,p1) + cross(r2,p2) + cross(r3,p3);
figure(5)
plot(L)
xlabel("timestep");
ylabel("angular momentum");
saveas(gcf, "unique_solutions/L_" + unique(i) + ".png");
end