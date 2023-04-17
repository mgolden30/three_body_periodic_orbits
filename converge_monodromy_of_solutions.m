%{
Let's see if we can converge the monodromy matrix in Fourier space
%}

clear

%load the particular solution you care about
load('solutions/3.mat');

addpath('functions');

T = y(end);
N = (numel(y)-1)/8;

y = reshape( y(1:8*N), [8,N] );

%subsample
%y = y(:,1:2:N);
%N = size(y,2);

Phi = zeros(8,8,N);

Phi(:,:,1) = eye(8);

I = eye(8);
h = T/N;

for i = 2:N
  J  = velocity_jacobian( y(:,i-1) );
  Jp = velocity_jacobian( y(:,i)   );

  Phi(:,:,i) = (I - h/2*Jp)\(( I + h/2*J ) * Phi(:,:,i-1));
end

%Take the final step
J  = velocity_jacobian( y(:,N) );
Jp = velocity_jacobian( y(:,1) );
Phi_f = (I - h/2*Jp)\(( I + h/2*J ) * Phi(:,:,N));



%% Let's fine-tune L


%Use this to estimate L
L = logm(Phi_f)/T;


phi = zeros(8,8,N);
for i = 1:N
  phi(:,:,i) = Phi(:,:,i) * expm( -L * h * (i-1) );
end

phi = reshape( phi, [64, N] );
%phi = [phi, phi];
plot( squeeze( real( phi(1,:) )) );




%% Newton-Krylov

z = zeros( 8*8*(N + 1), 1 );
z( 1:8*8*N ) = reshape( phi, [8*8*N,1] );
z( 8*8*N + (1:8*8) ) = L;

max_iterations = 32; %max number of Newton steps
every    = 4;
hookstep = 1;    %Fraction Newton step
restart  = 256;  %Krylov subspace size
tol      = 1e-9; %tolerance for GMRES
maxit    = 1;    %Only iterate Krylov subspace once

exit_condition = 1e-14; %Stop Newton at this error

normb = zeros( max_iterations, 1);

J = zeros(8,8,N);
for j = 1:N
  J(:,:,j) = velocity_jacobian( y(:,j) );
end

for i = 1:max_iterations
  tic

  b = monodromy_objective(z,y,J,T);
  normb(i) = norm(b);

  %Check for convergence
  if( normb(i) < exit_condition )
    normb = normb(1:i);
    fprintf('Converged...\n');
    break;
  end
  
  %finite difference the Jacobian
  h = 1e-6;
  A = @(v) ( monodromy_objective( z + h*v, y, J, T ) - b )/h;
  
  %Estimate Newton step with GMRES
  [step, ~] = gmres( A, b, restart, tol, maxit, [], [], 0*b );

  %Check that the NEwton step isn't too large
  hook = hookstep;
  step_size = norm(step)/norm(y);
  if step_size > hook
    hook = hook/step_size;
  end
  
  %If you are close to a solution, take a full Newton step
  if normb(i) < 1e-3
    hook = 1; 
  end
  
  %Actually take Newton step
  z = z - hook*step;
  
  %Print info to screen
  fprintf( 'step %d: |b| = %.3e, GMRES residual = %.3e\n', i, norm(b), norm( A(step) - b )/norm(b) );
  toc
end

%%

%load("1_mon.mat");

L2   = reshape( z(8*8*N + (1:64)), [8,8] );
phi2 = reshape( z(1:8*8*N), [8,8, N] );
phi  = reshape(phi, [8,8,N]);

plot( squeeze( abs(phi(1,4,:) - phi2(1,4,:))));
%plot( squeeze( real(phi2(1,1,:))));

return

%% Find Jordan form

[V,D] = eigs(L2, 8);
imagesc( abs(V'*V) )




v1 = f(y(:,1));
v1 = v1./norm(v1)

tol = 1e-9;
maxit = 1024;
v2 = lsqr( L2, v1, tol, maxit );

v3 = reshape( [0 1; -1 0] * reshape( y(:,1), [2,4] ), [8,1] );
v3 = v3/norm(v3);

v4 = lsqr( L2, v3, tol, maxit );

[V,D] = eigs(L2, 8);

U = V;
U(:, 5:8) = [v1,v2,v3,v4];
J = inv(U)*L2*U;

norm( V*D*inv(V) - L2 )
norm( U*J*inv(U) - L2 ) 


%% Look at spectra of phi

figure(4);
spectra = fft( reshape(phi2, [64,N] ), [], 2 )/N;
nice = abs(spectra);
xx = 1:(N/2-1);
nice = nice(:,xx);
scatter( xx, max(nice), 'filled' );
set( gca, 'YScale', 'log')
ylabel('max|{\bf \phi}_\omega|');
xlabel('\omega')
xlim([0, N/2-1]);
drawnow
set( gcf, 'Color' ,'w');
saveas(gcf, 'figures/converged_fourier_coeffs_phi.pdf');


