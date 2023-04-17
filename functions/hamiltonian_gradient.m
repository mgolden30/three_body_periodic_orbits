function dH = hamiltonian_gradient( x )
  %assume size(x) = [8,N]
  N = size(x,2);

  dH = zeros(8,N);
  r1 = x(1:2,:);
  r2 = x(3:4,:);
  p1 = x(5:6,:);
  p2 = x(7:8,:);
  
  r12 = r1-r2;
  r13 = 2*r1 + r2;
  r23 = 2*r2 + r1;

  d12 = vecnorm(r12);
  d13 = vecnorm(r13);
  d23 = vecnorm(r23);

  dH(5:6, :) = 2*p1 + p2;
  dH(7:8, :) = 2*p2 + p1;

  dH(1:2, :) =  r12./d12.^3 + 2*r13./d13.^3 +   r23./d23.^3;
  dH(3:4, :) = -r12./d12.^3 +   r13./d13.^3 + 2*r23./d23.^3;
end