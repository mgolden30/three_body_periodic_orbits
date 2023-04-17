function J = velocity_jacobian(x)
  %{
  PURPOSE:
  Compute the analytic Jacobian of a series of points in space x  
  %}
  
  J = zeros(8,8);
  
  %physics names
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
  


  J(1:4,5:8) = eye(4);
  
  I = eye(2);
  J(5:6,1:2) =  -  I./d12.^3 + 3*(r12*r12')./d12.^5 ...
                -2*I./d13.^3 + 6*(r13*r13')./d13.^5;
  J(5:6,3:4) =     I./d12.^3 - 3*(r12*r12')./d12.^5 ...
                -  I./d13.^3 + 3*(r13*r13')./d13.^5;
  
  J(7:8,1:2) =     I./d12.^3 - 3*(r12*r12')./d12.^5 ...
                -  I./d23.^3 + 3*(r23*r23')./d23.^5;
  J(7:8,3:4) =  -  I./d12.^3 + 3*(r12*r12')./d12.^5 ...
                -2*I./d23.^3 + 6*(r23*r23')./d23.^5;
end