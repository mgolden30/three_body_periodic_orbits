function [x] = converge_single_shooting( x0, h, N, max_iterations, verbose )
  %{
  PURPOSE:
  converge a periodic orbit from single-shooting in time. 

  INPUT:
  x0 - vector [r1;r2;p1;p2;T] of size [9,1]
  h  - finite difference parameter
  N  - number of timesteps
  max_iterations - max number of Newton steps
  verbose - true or false
  %}
  
  x = renormalize(x0);

  F = @(x) rk4step( x(1:8), x(9), N ) - x(1:8); %objective function
  
  J = zeros(10,9); %Jacobian and phase constraints
  Fx = zeros(10,1);
  for i = 1:max_iterations
    Fx(1:8) = F(x); %evaluate objective function

    %If verbose is true, print the current error
    if verbose
      fprintf( "Iteration %d: |F| = %.4e\n", i, norm(Fx) );
    end

    %Finite difference the Jacobian
    for j = 1:9
      x2 = x;
      x2(j) = x2(j) + h;
      J(1:8,j) = ( F(x2) - Fx(1:8) )/h;
    end
    J( 9,1:8) = f(x); %orthogonality to velocity
    J(10,1:8) = reshape( [0,1;-1,0]*reshape(x(1:8),[2,4]), [8,1] );

    %compute Newton step dx
    tolerance = 1e-4;
    [dx, ~] = lsqr(J,-Fx,tolerance);

    %make it smaller if needed
    dx = modify_newton_step(dx,x);
    
    %take the step
    x = x + dx;

    %fix energy H = -1
    x = renormalize(x);
  end
end


function dx = modify_newton_step(dx,x)
  % check to see if dx is too big
  % if it is too big, make is smaller

  r = 1e-2;
  if norm(dx) > r * norm(x)
    dx = dx/norm(dx)*norm(x)*r;
  end
end

function x = renormalize(x)
  %{
  Change normalization of energy
  %}
  
  H = hamiltonian(x);
  if ( H > 0 ) 
    %currently a scattering state, don't rescale anything. Hopefully we
    %will find our way back to negative energy.
    return;
  end

  l = sqrt(abs(H));

  x(1:4) = x(1:4)*l^2;
  x(5:8) = x(5:8)/l;
  x(9)   = x(9)*l^3;
end