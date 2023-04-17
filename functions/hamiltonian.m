function energy = hamiltonian( x )
  %Compute velocity of many points in phase space
  N = size(x,2);
  
  r1 = x(1:2,:);
  r2 = x(3:4,:);
  
  p1 = x(5:6,:);
  p2 = x(7:8,:);
  
  r12 = r1-r2;
  r13 = 2*r1 + r2;
  r23 = 2*r2 + r1;
  
  d12 = vecnorm( r12 );
  d13 = vecnorm( r13 );
  d23 = vecnorm( r23 );
  
  energy = sum( p1.^2 + p2.^2 + p1.*p2 ) - 1./d12 - 1./d13 - 1./d23;
end