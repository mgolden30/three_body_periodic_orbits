
addpath("functions/")

unique = [1,2,4,5,7,8,9,22,23,33,67,90,91,106,108,156];


for i = 1:numel(unique)
  load("solutions/" + unique(i) + ".mat");
  animate_state(y, "unique_solution_videos/" + unique(i) );
end