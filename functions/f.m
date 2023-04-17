function velocity = f( x )
  %Compute velocity of many points in phase space
  N = size(x,2);
    
  velocity = zeros( 8,N );
  
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
  
  velocity(1:2,:) = p1;
  velocity(3:4,:) = p2;
  
  velocity(5:6,:) = -r12./d12.^3 - r13./d13.^3;
  velocity(7:8,:) =  r12./d12.^3 - r23./d23.^3;

  %normalize velocity
  %velocity = velocity./vecnorm(velocity);
end