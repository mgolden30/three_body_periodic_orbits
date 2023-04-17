function Ju = jacobian_objective( x, u )
  uT= u(end);
  T = x(end);
  
  N = (numel(x)-1)/8;
  
  x = reshape( x(1:8*N), [8,N] );
  u = reshape( u(1:8*N), [8,N] );

  %compute velocity
  dvu = action_of_velocity_jacobian(x,u);
  v   = f(x);

  k = 0:(N-1);
  k( k> N/2) = k(k>N/2) - N;
  k(1) = -1i; %this way the zero mode is just multiplied by 1
  
  int_tau = @(v) real(ifft( fft(v,[],2)./(1i*k), [], 2));
%  v2 = fft(v, [], 2); %apply to second axis
%  v2 = v2./(1i*k); %integrate in time via Fourier 
%  int_v = real(ifft(v2, [], 2));

  int_v   = int_tau( v   );
  int_dvu = int_tau( dvu );
 
  
  Ju = zeros( [8*N+1,1] );
  
  u2 = u - mean(u,2); %subtract off zero mode
  
  %Take Jacobian analytically
  Ju(1:8*N) = reshape( u2 - int_dvu*T/(2*pi) - int_v*uT/(2*pi), [8*N,1] );
    
  %Add phase condition at the end to make Jacobian square
  dH = hamiltonian_gradient(x);
  Ju(end) = sum(sum(dH.*u))/N;
end