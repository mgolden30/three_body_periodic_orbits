function ob = monodromy_objective( z, y, J, T )
  N = (numel(z))/64 - 1;
  
  phi = reshape( z(1:8*8*N), [8,8,N] );
  L   = reshape( z( 8*8*N + (1:64)), [8,8] );

  v = zeros(8,8,N);
  for i = 1:N
    %J = velocity_jacobian( y(:,i) );

    v(:,:,i) = J(:,:,i)*phi(:,:,i) - phi(:,:,i)*L;
  end
  
  phi = reshape( phi, [8*8,N] );
  v   = reshape( v,   [8*8,N] );


  k = 0:(N-1);
  k( k> N/2) = k(k>N/2) - N;
  k(1) = -1i; %this way the zero mode is just multiplied by 1
  
  v2 = fft(v, [], 2); %apply to second axis
  v2 = v2./(1i*k); %integrate in time via Fourier 
  int_v = ifft(v2, [], 2); %Don't take the real part of this!!!!

  %int_v = real(int_v);

  % int_v is the inverse
  ob = zeros( [8*8*(N+1),1] );
  
  phi2 = phi - mean(phi,2); %subtract off zero mode
  
  ob(1:8*8*N) = reshape( phi2 - int_v*T/(2*pi), [8*8*N,1] );
  
  %Add phase condition at the end to make Jacobian square
  ob( 8*8*N + (1:64) ) = phi(1:64, 1) - reshape( eye(8), [64,1] );
end