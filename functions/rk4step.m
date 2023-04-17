function xp = rk4step( x, T, N )
  % Step forward in time by amount T for N steps
  % total time integrated will be dt*N

  xp = x;
  dt = T/N;
  for i = 1:N
    k1 = dt*f( xp        );
    k2 = dt*f( xp + k1/2 );
    k3 = dt*f( xp + k2/2 );
    k4 = dt*f( xp + k3   );
  
    xp = xp + (k1 + 2*k2 + 2*k3 + k4)/6;
  end
end