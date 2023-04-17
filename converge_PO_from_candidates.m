%{
Written by Matt Golden 2023

This is a variational solver for periodic orbits in the three-body problem
confined to a 2D plane.
%}

addpath('functions/');


clear;
clf

N = 1024;            % resolution in time
y = zeros( 8*N+1, 1 ); % trajectory vector
 
for a = 1:194

load( "candidates/" + a + ".mat" );

h = 1e-6;
N = 1024;
max_iterations = 32;
verbose = true;
[x] = converge_single_shooting( x, h, N, max_iterations, verbose )

y(1:8) = x(1:8);
y(end) = x(9);
y = generate_timeseries(y);
  
max_iterations = 1024; %max number of Newton steps
every    = 4;
hookstep = 0.05; %Fraction Newton step
restart  = 64;   %Krylov subspace size
tol      = 1e-9; %tolerance for GMRES
maxit    = 1;    %Only iterate Krylov subspace once
close_to_converged = 1e-2;

exit_condition = 1e-12; %Stop Newton when |b| is less than this
normb = zeros( max_iterations, 1);

for i = 1:max_iterations
  Hs = hamiltonian( reshape(y(1:8*N), [8,N] ) );
  if any( abs(Hs+1) > 1 )
    %failed
    break;
  end

  b = objective(y);
  normb(i) = norm(b);
  
  %Check for convergence
  if( normb(i) < exit_condition )
    normb = normb(1:i);
    fprintf('Converged...\n');
    save("solutions/" + a, "y");

    plot_state(y,1);
    drawnow;
    saveas(gcf, "solutions/figures/"+ a+ ".png");

    break;
  end
  
  %finite difference the Jacobian
  h = 1e-7;
  A = @(v) ( objective( y + h*v) - b )/h;
 

  %Estimate Newton step with GMRES
  [step, ~] = gmres( A, b, restart, tol, maxit );
  
  %Check that the NEwton step isn't too large
  hook = hookstep;
  step_size = norm(step)/norm(y);
  if step_size > hook
    hook = hook/step_size;
  end
  
  %If you are close to a solution, take a full Newton step
  if normb(i) < close_to_converged
    hook = 1;
  end
  
  %Actually take Newton step
  y = y - hook*step;
  
  %Print info to screen
  fprintf( 'step %d: |b| = %.3e, GMRES residual = %.3e\n', i, norm(b), norm( A(step) - b )/norm(b) );
  
  %Only make a plot every few steps
  if( mod(i-1,every) ~= 0 )
    continue; 
  end
  figure(1);
  plot_state(y,0);
  drawnow
end

end


function [d, Ts, ps] = search_for_orbits(N)
  N1 = 128; %points in momentum;
  N2 = 128; %points in T

  ps = linspace(-1, 1,N1);
  Ts = linspace( 4,13,N2);

  d = zeros(N1,N2);

  for i = 1:N1
    for j = 1:N2
      x = [1;0; 0;0; 0;0; ps(i);1; Ts(j)];
      x = renormalize(x);
      d(i,j) = norm( x(1:8) - rk4step(x(1:8), x(9), N) );
    end
  end

  imagesc(Ts, ps, d);
  xlabel("T_0");
  ylabel("p_0");
  axis square
 
  colorbar();
  caxis([0 2]);
end