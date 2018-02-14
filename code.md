---
layout: page
title: Code
permalink: /code/
---



Other than the ones with (in progress) in front of them, these are generally personal projects or very simplified toy models meant to explore a research idea and are not meant to be used for research. 
# Oceanography 

[**accesswoa13.m (in progress)**](https://github.com/FaizanHaque/matlab-tools/blob/master/accessWOA13.m)

Fetches the world ocean Atlas data using OpenDAP and if that's not available directly downloads the netCDF file.


# Agent/Individual based modelling 


**diversity.nlogo**

**brownian_bug.m**
I have a LxL box in which I initially generate N_0 "bugs" at uncorrelated x and y positions. Each of these N_0 bugs is assigned a color.
At every time step, each bug goes through a two step process:
1) The bug dies with probability q (it is multiplied by 0), the bug doubles with probability p (it is multiplied by 2) or the bug just survives with probability 1-q-p (multiplied by 1).
Any bug that reproduces has a duplicate at its exact same location and the duplicate is assigned the same color as the parent bug.
2) The x and y coordinates of each bugs is updated to simulate Brownian motion. This causes the colors (or the 'lineage') to spread.
The last part visualizes the spatial distribution of the bugs t=t0, t=t1 and t=tFinal where t0<t1<tFinal


**brownianGenetics.m**
Similar to above, except now instead of assigning colors, each bug is assigned a genome with Ng genes described by a 1 or 0. These are initially randomly assigned with probabilities and there is a set mutation per gene per time step.






