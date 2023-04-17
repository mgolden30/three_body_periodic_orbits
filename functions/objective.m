  function ob = objective( x )
  T = x(end);
  N = (numel(x)-1)/8;
  
  x = reshape( x(1:8*N), [8,N] );
  
  %compute velocity
  v = f(x);

  k = 0:(N-1);
  k( k> N/2) = k(k>N/2) - N;
  %k( k == N/2) = 0; %dealias
  k(1) = -1i; %this way the zero mode is just multiplied by 1
  
  v2 = fft(v, [], 2); %apply to second axis
  v2 = v2./(1i*k); %integrate in time via Fourier 
  int_v = real(ifft(v2, [], 2));

  % int_v is the inverse
  ob = zeros( [8*N+1,1] );
  
  x2 = x - mean(x,2); %subtract off zero mode
  
  ob(1:8*N) = reshape( x2 - int_v*T/(2*pi), [8*N,1] );
  
  %Add phase condition at the end to make Jacobian square
  %ob(end) = mean( (hamiltonian( x ) + 1).^2 );
  ob(end) = mean(hamiltonian(x)) + 1;
end