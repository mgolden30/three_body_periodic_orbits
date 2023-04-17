This code should allow one to find periodic orbits. To find these yourself, run the code in the following order.


1. hunt_for_candidates.m 
   This will perform random numerical integration to hunt for near-misses (George Carlin thinks they should be called near-hits)
   Results will be stored in the "candidates" directory

2. convergePO_from_candidates.m will do Newton on the single-shooting objective, and then Newton on the spectral representation. I have found this helps convergence quite a bit.
   Converged state will be stored in "solutions" with plots of trajectories in "solutions/figures"


