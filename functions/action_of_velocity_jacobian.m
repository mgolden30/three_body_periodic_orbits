function Ju = action_of_velocity_jacobian(x, u)
  %{
  PURPOSE:
  Compute the exact action of the velocity Jacobian
  %}

  N = size(x,2);
    
  Ju = zeros( 8,N );
  
  %physics names
  r1 = x(1:2,:);
  r2 = x(3:4,:);
  p1 = x(5:6,:);
  p2 = x(7:8,:);
  
  %physics names of perturbation vector
  ur1 = u(1:2,:);
  ur2 = u(3:4,:);
  up1 = u(5:6,:);
  up2 = u(7:8,:);
  
  r12 = r1-r2;
  r13 = 2*r1 + r2;
  r23 = 2*r2 + r1;
  
  ur12 = ur1-ur2;
  ur13 = 2*ur1 + ur2;
  ur23 = 2*ur2 + ur1;
  
  d12 = vecnorm( r12 );
  d13 = vecnorm( r13 );
  d23 = vecnorm( r23 );
  
  Ju(1:2,:) = up1;
  Ju(3:4,:) = up2;
  

  Ju(5:6,:) = -ur12./d12.^3 + 3*( sum(ur12.*r12) ).*r12./d12.^5 ...
              -ur13./d13.^3 + 3*( sum(ur13.*r13) ).*r13./d13.^5;
  Ju(7:8,:) =  ur12./d12.^3 - 3*( sum(ur12.*r12) ).*r12./d12.^5 ...
              -ur23./d23.^3 + 3*( sum(ur23.*r23) ).*r23./d23.^5;
end