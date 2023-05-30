# Running the code

For this worksheet, we have extended our previous work from the first worksheet, as to be able to simulate the energy flow (heat) in the liquid due to reservoirs being in contact with the fluid, and to run different geometry files, not only for the lid driven cavity (LidDrivenCavity) which has a very simple geometry: a box-shaped geometry where only the obstracles are in the edges of the box. Now the code can work for any arbitrary configuration of geometry, which can be provided by the user in a .pgm file. In order to make use of an arbitrary geometry, make sure to follow a simple set of rules:

- The file must have a squared array of integers from 0-8.
- Each of the digits represents a different type of tile. This convention for the assignment of different tiles must be followed for the simulation to work correctly (check the simulation section for more information on it).
- Specify all the required parameters as needed. If you want to make use of the heat equation make sure to initialize all the necessary parameters, the code will not check this for you.
- When creating arbitrary geometries, make sure that a boundary cell has not more than two (2) fluid neighboors. The code will let you know if you did something incorrectly, but will not terminate. We will implement try/catch blocks further down the semester as for the user to not have to constantly check the terminal in case he did something wrong so that he manually terminates the execution of the binary (or simply implement a method which checks the geometry before creating the simulation case and terminates if there is something wrong with it).

If you want to run different geometries, make sure the different files are in the example_cases folder. Each case must have the aforementioned .pgm file and a .dat file containing all the necessary parameters. If you want to create a new geometry, called for example mygeometry, make sure to name the two files the same way inside a folder. In such case, inside the example_cases there must be a directory called mygeometry containing two files: mygeometry.pgm and mygeometry.dat. Finally, for the program to recognize mygeometry as a valid geometry, edit in a text edit the main.cpp file. In the main, at the line nine (9) you will see a vector named cases, please make sure to add the name of your arbitrary geometry to it. In this case, you should add a new entrance to the vector with the string mygeometry.

In order to fork our repository and be able to produce results, make sure to run the following commands in your shell prompt:

shell
git clone https://gitlab.lrz.de/00000000014B328F/group-d-cfd-lab.git
cd group-d-cfd-lab
mkdir build
cd build
cmake ..
make


Afterwards, you should see a new binary file in your build directory named fluidchen. To execute the binary you have the option to run it for as many geometries as you want on the run, while just having to compile the code only once. In order to do that run the binary followed with the names of the geometries you want to run without any extension, separed by a blank space. If for example you want to run your newly defined geometry and the lid driven cavity, the way to run it would be with

shell
./fluidchen LidDrivenCavity mygeometry


This code will generate all the necessary output files for you to visualize the simulation in ParaView or any similar software of your preference. This process will create an output filder in ../example_cases/file_name/file_name_Output, which contains the .vtk files obtained from the simulation. If you want to use the heat equation, make sure that the variable energy_eq is set to on in the .dat file, otherwise, even if all the other parameters were correctly initialized, the simulation will not simulate the temperature.

# The Simulation

To correctly initialize the different tiles, follow the convention stated below:

0 : fluid.
1 : inflow.
2 : outflow.
3 : fixed wall (adiabatic).
4 : fixed wall (hot).
5 : fixed wall (cold).
8 : moving wall.

Note that tehcnically there is no difference between hot and cold walls, just their temperature which is input by the user in the .dat file anyway, but if possible try to maintain this convention to avoid possible confussion. The extra id's 6 and 7 can be used by the user to define extra custom walls.

To simulate the energy transfer in the fluid, other than using the Navier-Stokes equations to simulate the dynamics of the fluid, we have used the heat equation which is given by:

math
\dfrac{\partial T}{\partial t} + \vec{u}\dot\vec{\nabla}T = \alpha\Delta T + Q,

where T denotes the temperature, time is t, $\vec{u}$ is the velocity vector, $\vec{\nabla}$ is the gradient, or vector of derivatives, Q is the extenal sources/sinks of energy and $\Delta$ is the laplacian (summation over the second derivatives).

Additionally, we have used the Boussinesq approximation, for which we assume that the density remains constant through the whole fluid (we are assuming incompressible fluids) except for terms which include buoyancy, which are terms coming from the fact that materials tend to expant (their density decreases) when they heat up.

Finally, to ensure the different observables are held constant as expected from the input variables given by the user, we have used the following boundary conditions for the different type of cells:

- Hot/cold walls: Since they need to have a constant temperature throughout the whole simulation, we have used Diriclet boundary conditions as to ensure they always retain their temperature T.
- Abiadatic walls: As the name implies, for this type of wall, we have used Neumann conditions with zero derivative, meaning adiavatic conditions, which means that there is no grandient of temperature between said wall and its neighboor, meaning that both of them have the same temperature.
- Inflow: For this boundary tile we used the fact that naturally, their velocity must be constant, as we used Diriclet boundary conditions for their velocity, using the input value given by the user, and for the pressure we again used adiabatic conditions as their pressure should not affect that of the neighbooring cells.
- Outflow: Finally, for the outflow cells we used Dirichlet boundary conditions for the pressure, to maintain it as zero. Although this value should be set to one atmosphere, which is the pressure assumed outside the fluid container, this does not affect either the simulation or the visualization results, as we are more concern about the movement of the fluid itself than how it flows out. For the velocity we have used Neumann conditions with no gradient again, since naturally, is the fluid is moving at high speeds just beside the opening, then we would physically expect the fluid to come out at high speed.

# Challenges

Finally, on top of having to think of an efficient way of reading and implementing the new changes to the workflow, we found it particularly difficult to debug the code, as we had some small bugs passing over from the first worksheet, which were very subtle and difficult to notice, and took us some time trying to test and debug the code as for it to work. This was reflected on the residual of the SOR method, which actually had nothing to do with the actual bug, as we did not had to change the implementation of the SOR at all, which also misguided us into what the bug might have been. Finally, implement also the boundary conditions of the inflow and outflow cells was a complicated at the beginning, as there were no instructions provided on out to deal with them, so we had to use our own knowledge and physical intuition to do something as for the simulation to be able to fully and correctly run, so that we could visualize the results and decide on the go what to change.

# Results and Discussion
## 1.3 and 1.4

All the tasks mentioned in `Section 1.3 and 1.4` of the worksheet have been implemented and the visualizations of the example cases are shown below for reference. Using the base parameters from 1.4, the converged solution was visualized and verified against the posted solution on Moodle. The results are similar to the sample solution, thus verifying the correctness of our code.

## The Lid Driven Cavity
<div align="center">
  <img width="800" height="550" src="Cavity_Velocity.png">
  <figcaption>Velocity Surface Plot</figcaption>
</div>
One of the requirements of this worksheet was to preserve what was already implemented previously. As can be seen from the velocity surface plot for lid driven cavity, the accuracy of results in this case is maintained.

## Task 1.4(b) The Karman Vortex Street 
<div align="center">
  <img width="800" height="550" src="ObsFlowVel.png">
  <figcaption>Velocity Surface Plot</figcaption>
</div>
The fluid inflow was set to u = 1.0 and v = 0.0. The upper and Lower boundaries have no slip boundary conditions imposed. The results from the surface plot is consistent with the sample output. The velocity has a minimum value of 0 at the boundaries and also right behind the areas of the obstacle.

## Task 1.4(c) Channel flow with a backward facing step 
<div align="center">
  <img width="800" height="550" src="BFS_Velocity.png">
  <figcaption>Velocity Surface Plot</figcaption>
</div>
The inflow conditions in this case was set to u = 1.0 and v = 0.0. The boundary conditions imposed where no-slip conditions at the upper and lower walls. The geometry file creates an onstacle domain representing a square filling up half of the channel height. The results as can be seen from the velocity surface plot seems to be consistent with the sample output.

## Task 2.2(d) Natural Convection
<div align="center">
  <img width="800" height="550" src="ConvectionVel.png">
  <figcaption>Velocity Surface Plot</figcaption>
</div>
In this case, the heat equation was used to calculate the temperature values. As can be seen from the velocity surface plots, the fluid appears to have a circular motion in between the hot and cold walls. It is also interesting that the velocity is zero around the center of the domain.

<div align="center">
  <img width="800" height="550" src="ConvectionTemps.png">
  <figcaption>Temperature Surface Plot</figcaption>
</div>
The temperatures at the boundaries where applied as boundary conditions. The results shown here where obtained with the first set of values provided in the worksheet with a dt value of 10. As expected, the temperature values at the end of the simulation shows that the heat flow is infact from the hot walls towards the cold and insulated walls.

## 2.1 and 2.2

All the tasks mentioned in `Section 2.1 and 2.2` of the worksheet have been implemented. However, some issues were encountered during the simulation. The residulas are diverging giving a value of Inf or NaN for RayleighBenard case and FluidTrap case. 

# Challenges

1) Ensuring that the Boundary conditions are correctly declared and defined for new Boundary classes.

2) Understandng the cause of various errors and rectifying them.

3) Implementing the energy equation on/off condition for multiple cases in thr files. 





