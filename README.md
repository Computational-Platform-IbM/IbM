# IbM

Individual based model on bacterial colonies

Mathemathical model simulating the diffusion of compounds and bacterial growth.

To execute the model (Granule versions):

1. Create .mat files using 'mat_creator(number_replicates)'
2. Drag the desire sim\_####.mat to main folder (where IbM.m is found)
3. Execute the model with a call to 'IbM(sim_number)'
4. To visualise the bacteria of the results, use 'python plotBacsFast.py (sim_number) (-nf)' in command prompt (Anaconda prompt or others).
   Add -nf only if the simulation has not finished yet.

More information on the background of the model can be found in the Materials & Methods file
