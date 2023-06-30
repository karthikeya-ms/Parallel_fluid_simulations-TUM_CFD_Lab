## Worksheet 2

### Plain Shear Flow
We begin with modeling the simple case of Plain Shear Flow in a pipe. No slip boundary conditions are used in this case. It can be observed that the pressure is continuosly decreasing along the length of the pipe. The reason for the decrease in pressure is frictional loss due to the presence of viscosity. From velocity contour, it can be observed that the highest velocity is seen along the centerline after the boundary layer is fully developed and a steady flow is reached. Velocity at wall is zero as expected. The velocity of 1.5 in the center matches with the value obtained from the analytical expression of velocity profile in pipe flow. 

<p>
<img src= "docs/Worksheet2_Plots/PlainShear_Pressure_Contour.png" width="500"> | <img src="docs/Worksheet2_Plots/PlainShear_VelocityX_Contour.png" width="500">
</p>
<em>Pressure Contour plot for pipe flow [left], and Velocity contour plot [right] </em>

<p>
<img src= "docs/Worksheet2_Plots/PlainShear_VelocityX_Profile.png" width="500"> | <img src="docs/Worksheet2_Plots/PlainShear_Pressure_Profile.png" width="500">
</p>
<em>Velocity profile along the section [left], and Pressure and Velocity magnitude profile along the centerline [right] </em>


### Karman Vortex Simulation
Karman vortex refers to vortex shedding caused by flow around a blunt obstacle. In present simulation, the fluid enters from left with a constant velocity and encounters a tilted plate. This leads to vortex shedding which can be seen in the velocity or pressure contour. A small recirculation zone behind the obstacle is observed as expected. You can find link to the animation file [here](https://gitlab.lrz.de/00000000014ADC3D/cfd-lab-group-mib/-/blob/Worksheet2/docs/Worksheet2_Plots/KarmanVortex_Velocity.avi)

<p>
<img src= "docs/Worksheet2_Plots/Obstcle_VelocityMag_Contour.png" width="500"> 
</p>
<em>Velocity Contour Plot for flow around the obstacle</em>

### Flow Over a Step
In this case, we simulate the effect of a sudden step in the flow path. The fluid enters from left with a constant velocity. No slip boundary conditions are used in the upper and lower plates. The flow expands on encountering the step. A recirculation zone is formed just downstream of the step. The velocity is higher near the entrance owing to the smaller area to conserve the mass flow rate. However, the velocity is more of less unifrom across the section. In the region towards the outlet, velocity is lower in magnitude and non-uniform across the section (parabolic as seen in the pipe flow case above).

<p>
<img src= "docs/Worksheet2_Plots/Step_Pressure_Contour.png" width="500"> | <img src="docs/Worksheet2_Plots/Step_VelocityMag_Contour.png" width="500">
</p>
<em>Pressure Contour Plot [left], Velocity Contour Plot [right]</em>


### Natural Convection
This is the first case where we simulated a flow sololey due to gravity without prescribing any fixed velocity at the boundary. The simulation domain is a square where top and bottom walls are insulated, left wall is heated and right wall is cooled. Gravity acts on the system.
We conducted a first simulation with viscosity nu = 0.001 and thermal diffusivity alpha = 0.000142857. In the following picture we can notice how the hot fluid tends to expand more in the higher region of the domain, while the cold fluid expands more in the lower region. Hot fluid is less dense and hence it tends to rise, while on the other hand cold fluid is denser and tends to descend. This is exactly what we can observe through the velocity contour. Glyphs on the left-hot side point upwards, while those on the right-cold side point downwards.

<p> 
<img src="docs/Worksheet2_Plots/NaturalConvection_Temp_contour_velGlyph_case1.png" width="500" >
</p>
<em>Temperature Contour plot superimposed with Velocity Glyphs, nu = 0.001 and alpha = 0.000142857.</em>

To understand why the flow moves with a clockwise circular motion we also need to have a look at the pressure profile. Pressure is overall higher in the top region of the domain, and lower on the bottom part. This is due to buoyancy effect. Since hot gas moves upward relative to cold, the region with more hot gas (top) will see pressure building up at top. For both top and bottom region we can still see a horizontal difference in pressure. On the top region pressure is higher on the left corner, while in the bottom region pressure is higher on the right corner. Fluid moves from regions of higher pressure to those of lower pressure, hence in the top region of the domain, gas moves from left to right and in the bottom region it moves from right to left. 

<p> 
<img src="docs/Worksheet2_Plots/NaturalConvection_Press_contour_velGlyph_case1.png" width="500">
</p>
<em>Pressure Contour plot superimposed with Velocity Glyphs, nu = 0.001 and alpha = 0.000142857. </em>

The overall effect is a circular motion. Motion along the y-axis is temperature driven, while motion along the x-axis is pressure driven. Pressure profile is given by the temperature we prescribed at the boundaries.

The same observations can be also made for the second simulation, which we conducted at nu = 0.0002 and alpha = 0.000028571. However, we can point out some substantial differences.
Due to the lower thermal diffusivity in the second simulation, heat is conducted less quickly and hence heated gas moves less rapidly. This is well showed in the temperature profile, here the region of heated fluid is much more schrinked if compared to simulation one. This also reflects in the pressure profile, where we can observe an over all lower value of pressure. Pressure has a 0 value on most of the domain and presents positive values only in the top left corner.The lower value in viscosity results in fluid layers moving more easily on top of each other. 
As an overall result, the vortex center has now moved to the left side of the domain, while in simulation one it was at the exact center of the domain.

<p> 
<img src="docs/Worksheet2_Plots/NaturalConvection_Temp_contour_velGlyph_case2.png" width="500"> | <img src="docs/Worksheet2_Plots/NaturalConvection_Press_contour_velGlyph_case2.png" width="500">
</p>
<em>Temperature Contour plot superimposed with Velocity Glyphs [left] and Pressure Contour plot superimposed with Velocity Glyphs [right], nu = 0.0002 and alpha = 0.000028571. </em>


### Fluid Trap Simulation
In this case we simulate a heat driven flow with obstacles. You can observe from the Temperature contour plot that the left wall is the hot wall while the right wall is the cold wall. The flow reaches an equilibrium where the hot fluid from left section is not able to propagate into the cold section on right and cold gas from right section is not able to reach the left section. Further, in the pressure contour, one may not the high pressure region in top-left and bottom-right of the domain. This is because the hot gas in left section moves upwards creating a high pressure region top-left of domain, while high amount of cold gas sinks to bottom in the right section creating high pressure zone. If we observe the velocity glyphs, we can see that some amount of hot gas from left section reaches right and moves upward, causing slight higher temperature in top-right compared to bottom-right. Similar observation can be made for left section. 

<p> 
<img src="docs/Worksheet2_Plots/FluidTrap_TempContour_VelGlyph.png" width="500">
</p>
<em>Temperature Contour plot superimposed with Velocity Glyphs</em>

<p> 
<img src="docs/Worksheet2_Plots/FluidTrap_PressContour_VelGlyph.png" width="500">
</p>
<em>Pressure Contour plot superimposed with Velocity Glyphs</em>

### Rayleigh Benard Convection

Rayleigh Benard is a natural convection occuring in a horizontal fluid heated from bottom. The fluid in contact with the hot wall at bottom rises up till it reaches the top (typically a fluid-fluid interface). The top boundary is modelled as a cold wall in present simulation. The fluid at top turns to side but cannot move far as it encounters another coloum of fluid around it. It instead forms a loop with itself creating isolated cells called Benard Cells. This is illustrated in the figure below.

<p> 
<img src="docs/Worksheet2_Plots/ConvectionCells.png" width="500">
</p>
<em>Convection Cells or Benard Cells (Source: Wikipedia)</em>

In both the temperature and velocity contour plot, we can observe the formation of convection cells. Glyphs highlight the closed velocity loops as expected for the problem. However, the pattern is not very regular in the present case. This may be due to coarse grids or smaller domain size as Rayleigh Benard convection requires the height of the fluid layer to be much smaller than the horizontal dimensions. The simulation was repeated with a longer domanin but results obtained were qualitatively similar.

<p> 
<img src="docs/Worksheet2_Plots/RB_TempContour_VelGlyph.png" width="500">
</p>
<em>Temperature Contour plot superimposed with Velocity Glyphs</em>

<p> 
<img src="docs/Worksheet2_Plots/RB_VelContourGlyph.png" width="500">
</p>
<em>Velocity Magnitude Contour and Glyphs</em>

<p> 
<img src="docs/Worksheet2_Plots/RB_VelContourGlyph_longerdomain.png" width="500">
</p>
<em>Velocity Magnitude Contour and Glyphs for case with longer domain</em>

---

## Worksheet 1
### Pressure Visualization

The image below shows the pressure contour for the lid driven cavity simulation, after some thousands time steps. As the picture shows, the highest pressure is in the upper right corner of the lid driven cavity while the lowest pressure is observed at the left upper corner. If one also takes a look at the glyphs in the velocity profile (next section), then it's obvious that the fluid moves from high pressure to low pressure. The only exception is the top of the cavity. Here the cells used for the discretization are "attached" to the lid and hence move with it (from low pressure to high pressure). 

<p> 
<img src="example_cases/LidDrivenCavity/Plots/contour_pressure_gray.png" width="500"> | <img src="example_cases/LidDrivenCavity/Plots/contour_pressure_gray_scaled.png" width="500">
</p>
<em>Pressure Contour default [left] and scaled [right]</em>

### Velocity Visualization

A expected, velocity is zero at left, right and bottom boundaries, while it's equal to the lid's velocity (one) at the top boundary, where cells are attached to the lid (as explained in the previous section).

<p> 
<img src="example_cases/LidDrivenCavity/Plots/contour_u_gray.png" width="500">
<img src="example_cases/LidDrivenCavity/Plots/contour_v_gray.png" width="500">
<img src="example_cases/LidDrivenCavity/Plots/Glyph_velocity_gray.png" width="500">
</p>
<em>Velocity Contour in x-direction [left], Velocity Contour in y-direction [middle] and Velocity Glyph [right]</em>


### Examination of SOR solvers behavior depending on ω

To examine the SOR solver's behaviour depending on the relaxation factor ω, we conducted the simulation with different values of ω and for each one of these and for each time step, we searched for the maximum number of iterations necessary for the solver to converge (res < eps). Among the ω values we used, the one that required less iterations to converge was ω=1.9. For very low values of ω a very large number of iterations is necessary to obtain convergence. That's consistent with the SOR solver formulation since p(n+1)->p(n) for ω->0. In the table below we report also the time step at which the max iteration is reached, i.e. the time step that required the maximum number of iterations to converge.

ω | max-iter | timestep |
--- | --- | --- |
0.5 | 3499 | 592 |
1.0 | 1132 | 1992 |
1.3 | 655 | 1992 |
1.6 | 346 | 1992 |
1.7 | 260 | 1992 |
1.8 | 176 | 1992 |
1.9 | 154 | 1992 |
1.95 | 295 | - |
1.99 | 1350 | - |


### The algorithm’s behavior depending on Timestep (δt)

In our code we provided an implementation for adapting the time step size δt in accordance to the stability condition in Equation (13). However, in order to analyse the algorithm's behaviour depending on δt, we kept it constant and equal to the values reported in the table below throughout the whole simulation. We observed that the solution diverged for the timesteps 0.01, 0.03 and 0.05 and converged for timesteps 0.005, 0.007 and 0.009.

dt | timestep | stabilitiy |
--- | --- | --- |
0.05 | 5 | Diverged |
0.03 | 7 | Diverged |
0.01 | 71 | Diverged |
0.009 | - | Converged |
0.007 | - | Converged |
0.005 | - | Converged |

### The algorithm’s behavior depending on i_max = j_max

We also studied the algorithm's behaviour when using different values of i_max an j_max. The bigger i_max=j_max gets, the finer the grid gets (smaller dx and dy). To ensure the stability condition (Equation 12), δt also has to smaller. We conducted our simulations for varying values of i_max=j_max, but always keeping δt constant to 0.05. We then looked at u(i_max/2, 7*j_max/8) at the end of the simulation for each one of these time steps. The only converged result we obtained was for i_max=j_max=16 for kinematic viscosity, nu=0.01, in all the other cases, solution diverged. However, for nu = 0.001, we obtained convergence upto i_max=j_max=32. The solution diverged for cases with further refinement. In order to obtain physical results also for bigger values of i_max=j_max we would need to decrease the time step size δt. 

imax/jmax | Cell(i,j) | u(i,j) for nu=0.01 | u(i,j) for nu=0.001 |
--- | --- | --- | --- |
16 | (8,14) | 0.183833 | 0.139114 |
32 | (16,28) | -nan | 0.208146 |
64 | (32,56) | -nan | -nan |
128 | (64,112) | +nan | -nan |
256 | (128,224) | +nan | +nan |



### Effect of kinematic viscosity

We investigated the effect of kinematic viscosity on our simulation. We looked at the velocity profile for different values of kinematic viscosity (nu). For lower values of viscosity the effect of the moving lid can propagate deeper into the fluid, however the velocity magnitude is very high only for regions of fluid that are very close to the lid. This can be understood by considering the Reynolds number (Re), which is inversely proportional to the kinematic viscosity. At large Re, the boundary layer attached to the lid is similar to a turbulent boundary layer (zero relative velocity close the wall which decreases sharply as you move away from it). Also, large Re corresponds to greater mixing in the flow causing the effect of moving lid to be felt at distant locations in the cavity.

For large viscosity (low Re), the region of "high" (close to 1) velocity magnitude extends more (as in a laminar flow boundary layer) but fluid away from the lid remains largely unaffected (the rest of the cavity is at values very close to zero).

nu = 0.01 :
<img src="/example_cases/LidDrivenCavity/Plots/NuComparison/ustream_100_nu01.png" width="500"> | <img src="/example_cases/LidDrivenCavity/Plots/NuComparison/u_100_nu01.png" width="500"> 

nu = 0.002 :
<img src="/example_cases/LidDrivenCavity/Plots/NuComparison/ustream_100_nu002.png" width="500"> | <img src="/example_cases/LidDrivenCavity/Plots/NuComparison/u_100_nu002.png" width="500"> 

nu = 0.0005 :
<img src="/example_cases/LidDrivenCavity/Plots/NuComparison/ustream_100_nu0005.png" width="500"> | <img src="/example_cases/LidDrivenCavity/Plots/NuComparison/u_100_nu0005.png" width="500"> 

nu = 0.0001 :
<img src="/example_cases/LidDrivenCavity/Plots/NuComparison/ustream_100_nu0001.png" width="500"> | <img src="/example_cases/LidDrivenCavity/Plots/NuComparison/u_100_nu0001.png" width="500"> 
