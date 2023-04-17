function plot_state(x, closed)
  if size(x,2) == 1
    N = ( numel(x)-1 )/8;
    x = reshape( x(1:8*N), [8,N]);
  end

  r1 = x(1:2,:);
  r2 = x(3:4,:);
  r3 = -r1-r2;
  
  if closed == 1
    p = @(x,c) plot( [x(1,:), x(1,1)], [x(2,:), x(2,1)], 'LineWidth', 3, 'Color', c );
  else
    p = @(x,c) plot( x(1,:), x(2,:), 'LineWidth', 3, 'Color', c );
  end
  p(r1, 'red');
  hold on
    p(r2, 'blue');
    green = [8 125 10]/256;
    p(r3, green);
  hold off
  
  xlim([-1 1] * 3);
  ylim([-1 1] * 3);
  axis square
end