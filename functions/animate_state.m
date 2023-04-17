function animate_state(x, output)
  if size(x,2) == 1
    N = ( numel(x)-1 )/8;
    T = x(end);
    x = reshape( x(1:8*N), [8,N]);
  end

  %Period double the solution for the video
  x = [x,x];
  T = 2*T;
  N = 2*N;

  p = @(x,c) plot( x(1,:), x(2,:), 'LineWidth', 3, 'Color', c );
  
  r1 = x(1:2,:);
  r2 = x(3:4,:);
  r3 = -r1-r2;
  
  %fix rotation?
  v = sum(r1.^2).*r1 + sum(r2.^2).*r2 + sum(r3.^2).*r3;
  v = mean(v,2);
  norm(v)

  if norm(v) < 1e-7
    %the method with v will not work. Try a higher order tensor
    q = 2;
    v11 = sum(r1.^q).*r1(1,:).*r1(1,:) + sum(r2.^q).*r2(1,:).*r2(1,:) + sum(r3.^q).*r3(1,:).*r3(1,:);
    v12 = sum(r1.^q).*r1(1,:).*r1(2,:) + sum(r2.^q).*r2(1,:).*r2(2,:) + sum(r3.^q).*r3(1,:).*r3(2,:);
    v22 = sum(r1.^q).*r1(2,:).*r1(2,:) + sum(r2.^q).*r2(2,:).*r2(2,:) + sum(r3.^q).*r3(2,:).*r3(2,:);

    V = [mean(v11,2), mean(v12,2);
         mean(v12,2), mean(v22,2)];

    [U,D] = eigs(V);
    v = U(:,1);
    V
  end
  R = [v(1), v(2); -v(2), v(1)]/norm(v);

  r1 = R*r1;
  r2 = R*r2;
  r3 = R*r3;

  p = @(x,c) plot( x(1,:), x(2,:), 'LineWidth', 3, 'Color', c );
  
  vidObj = VideoWriter(output);
  vidObj.FrameRate = 30;
  open(vidObj);

  %I want this in real time, so figure out step n to get 30 fps
  steps = 30*T;
  for i = 1:steps
    e = round( i*N/steps ); %end timestep for this frame

    p(r1(:,1:e), 'red');
    hold on
    p(r2(:,1:e), 'blue');
    green = [8 125 10]/256;
    p(r3(:,1:e), green);
    
    %add dots for particles
    scatter( r1(1,e), r1(2,e), 100, 'filled', 'MarkerFaceColor', 'black' );
    scatter( r2(1,e), r2(2,e), 100, 'filled', 'MarkerFaceColor', 'black' );
    scatter( r3(1,e), r3(2,e), 100, 'filled', 'MarkerFaceColor', 'black' );
    hold off

    xlim([-1 1] * 3);
    ylim([-1 1] * 3);
    xticks([]);
    yticks([]);
    set(gcf, 'Color', 'w');
    set(gca, 'Visible', 'off');
    axis square

    drawnow
    
    % Write each frame to the file.
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  
  % Close the file.
  close(vidObj);
end