# dynamic-brainstates
This repository contains code to cluster BOLD timeseries data into clusters and to calculate brain dynamics parameters (FO = fractional occupancy, DT = dwell time, AR = appearance rate), based on the paper by Cornblath et al.:
Cornblath, Eli J., Arian Ashourvan, Jason Z. Kim, Richard F. Betzel, Rastko Ciric, Azeez Adebimpe, Graham L. Baum, Xiaosong He, Kosha Ruparel, Tyler M. Moore, Ruben C. Gur, Raquel E. Gur, Russell T. Shinohara, David R. Roalf, Theodore D. Satterthwaite, and Danielle S. Bassett (2020), “Temporal sequences of brain activity at rest are constrained by white matter structure and modulated by cognitive demands,” Communications biology, 3 (1), 261.


The code to cluster BOLD timeseries, calculate centroids, and generate radial plots is from the Cornblath Repo: https://github.com/ejcorn/brain_states
The code to generate brain maps of clusters is on Keith Jamison's repo: https://github.com/kjamison/atlasblobs



processing and denoising:
- colossus: /mnt/shared_data3/emo4002/preprocessing_conn/sf_processing.m and /f_denoising.m
- allsubs/ is where denoised functional output can be found
- functional connectivity extraction: extract_voxelwise_ts.m

