# Running the code

The goal of this worksheet is to implement parallelization of the already implemented previous worksheet using the MPI message passing library and make the simulation faster. To run the code use the following command - 

- mpirun -np _no._of_processors ./fluidchen Simulation_case_name

For eg: mpirun -np 4 ./fluidchen ChannelWithBFS

Doing this for the first time will prompt the program to work in sequential mode and will ask the user if they are sure to continue running it in sequential mode. Answer 'yes' or 'no'. If no, enter the number of processors in x and y direction and the program will save that information in the DAT file. Now, run the above command with the desired size again like so -

mpirun -np 4 ./fluidchen ChannelWithBFS

 Now the program should begin the desired simulation with the desired number of processor cores.

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

```
\dfrac{\partial T}{\partial t} + \vec{u}\dot\vec{\nabla}T = \alpha\Delta T + Q,
```

where T denotes the temperature, time is t, $\vec{u}$ is the velocity vector, $\vec{\nabla}$ is the gradient, or vector of derivatives, Q is the extenal sources/sinks of energy and $\Delta$ is the laplacian (summation over the second derivatives).

Additionally, we have used the Boussinesq approximation, for which we assume that the density remains constant through the whole fluid (we are assuming incompressible fluids) except for terms which include buoyancy, which are terms coming from the fact that materials tend to expant (their density decreases) when they heat up.

Finally, to ensure the different observables are held constant as expected from the input variables given by the user, we have used the following boundary conditions for the different type of cells:

- Hot/cold walls: Since they need to have a constant temperature throughout the whole simulation, we have used Diriclet boundary conditions as to ensure they always retain their temperature T.
- Abiadatic walls: As the name implies, for this type of wall, we have used Neumann conditions with zero derivative, meaning adiavatic conditions, which means that there is no grandient of temperature between said wall and its neighboor, meaning that both of them have the same temperature.
- Inflow: For this boundary tile we used the fact that naturally, their velocity must be constant, as we used Diriclet boundary conditions for their velocity, using the input value given by the user, and for the pressure we again used adiabatic conditions as their pressure should not affect that of the neighbooring cells.
- Outflow: Finally, for the outflow cells we used Dirichlet boundary conditions for the pressure, to maintain it as zero. Although this value should be set to one atmosphere, which is the pressure assumed outside the fluid container, this does not affect either the simulation or the visualization results, as we are more concern about the movement of the fluid itself than how it flows out. For the velocity we have used Neumann conditions with no gradient again, since naturally, is the fluid is moving at high speeds just beside the opening, then we would physically expect the fluid to come out at high speed.

# Challenges

The challenges with regards to this worksheet was working with MPI. implementing the communication methods and making sure that no deadlocks occur was very difficult. We also had to work with multiple sub domains and intitialize all the indeced and sizes very carefully. Figuring out how we communicate between the neighbouring threads and finding which threads to communicate to was also very challenging. Another difficulty was trying to figure out how to parse the geometry file and then decompose the geometry data and assign it to each seperate threads. Finally, assigning the local sub domain sizes was also really difficult. It was necessary to figure out the relationship between `iproc`, `jproc` and the subdomain sizes and assign them appropriately. 

# Results and Discussion
## Issues with the code
As can be seen after running our code, the intended result is not achieved. Here is a breakdown of the tasks and our approach at tackling them - 

### Tasks 1 and 2

In `case.cpp` we have implementations for reading the additional communication specific parameters. We also make sure that the size of the parallel process is the product of `iproc` and `jproc` and if this is not the case, we throw an error message regarding this. The program then goes on to decompose the domain. We implemented if-else conditions to assign the local max indeces of the sub domain. We check if `iproc` and `jproc` are multiples of `imax` and `jmax` (all processors take same number of tiles in x and y directions). If this is not the case, we create seperate cinditions for the appropriate assignment. After assigning the local max indeces, the doamin struct is initialized and the build_domain() method is called to build the sub-domain for each process. 
In Grid.cpp, the master rank reads through the whole geometry file. It then sends the geometry data in sizes of the subdomain to all other processes. The other processes recieve the geometry data in sizes of their subdomain and stores them in their local geometry data vector. The domain decomposition and assignment of the geometry data to the various sub domains is now complete.

### Task 3

In `communication.cpp`, we create a method, void communicate() which is the crux of all the ghost cell exchange between the threads. The method calculates the ranks of its right, left, top and bottom neighbour threads. If any of them is a solid boundary, then they are assigned MPI_PROC_NULL. Now, we determine the positioning of the sub domain relative to its neighbours. Based on this information, we initiate the appropriate communication methods. For instance, if a subdomain is at the left bottom corner, it only has a neighbour on the top and right sides. So it only sends and recieves from the right side. The method just accepts a single matrix as argument. MPI_Sendrecv is used for all the exchanges to avoid potential deadlocks. 

### Tasks 4, 5 and 6

In `cell.cpp` we create a boolean attribute and method as a flag to see if a particular cell is part of the ghost layer. In `grid.cpp`, we assign this boolean to true if it is indeed a ghost cell and we do not push them to the vector of fluid cells.
In `Fields.cpp` and `PressureSolver.cpp` we calculate the fluxes, temperature, pressure and velocities for the invidual subdomain and then send the resulting matrix to the communicate method in `communication.cpp` where the ghost layer cell data is communicated. In `case.cpp`, we calculate the local maximum velocities. The master rank accumulates all the local max velocities and determines the global max velocities. This value is the broadcasted to all the other processes. The global max velocity is used to determine the value of dt. 

### Tasks 6, 7 and 8

The output methods where modified as per the requirements of the worksheet.

### The Issues

The code does work for single processor. However, it does not work when used in parallel mode. This implies that the issue must be arrising from the communication methods. However we are unable to determine the exact nature of the issue. 






