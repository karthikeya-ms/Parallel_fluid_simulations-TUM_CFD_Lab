#include "Case.hpp"
#include "Enums.hpp"
#include "Discretization.hpp"

#include <algorithm>
#ifdef GCC_VERSION_9_OR_HIGHER
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <mpi.h>
#include <cstdlib>

#ifdef GCC_VERSION_9_OR_HIGHER
namespace filesystem = std::filesystem;
#else
namespace filesystem = std::experimental::filesystem;
#endif

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>
#include <limits>

Case::Case(std::string file_name, int argn, char **args) {
//Case::Case(std::string file_name) {
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);
    double nu;      /* viscosity   */
    double UI;      /* velocity x-direction */
    double VI;      /* velocity y-direction */
    double PI;      /* pressure */
    double GX;      /* gravitation x-direction */
    double GY;      /* gravitation y-direction */
    double xlength; /* length of the domain x-dir.*/
    double ylength; /* length of the domain y-dir.*/
    double dt;      /* time step */
    int imax;       /* number of cells x-direction*/
    int jmax;       /* number of cells y-direction*/
    double gamma;   /* uppwind differencing factor*/
    double omg;     /* relaxation factor */
    double tau;     /* safety factor for time step*/
    int itermax;    /* max. number of iterations for pressure per time step */
    double eps;     /* accuracy bound for pressure*/
    double UIN{0};     /* x-dir. velocity of fluid coming into the domain */
    double VIN{0};     /* y-dir. velocity of fluid coming into the domain */
    std::string energy_eq; /* boolean determining if the heat equation for computing temperature is to be used */
    double TI{0};      /* initial temperature */
    double beta{0};    /* thermal expansion coefficient */
    double alpha{0};   /* thermal diffusivity */
    int num_walls{6};  /* number of different wall classes (wall id ranging from 3-7). Initialized with arbitrary value to allocate heap memory for temps variable */
    double temps[num_walls]; /* Vector of temperatures of the different walls */
    std::string geo_file{"NONE"};     /* String with the name of the geometry file loaded */
    int iproc{1};    /* number of processors used to parallelize simulation in x-axis */
    int jproc{1};    /* number of processors used to parallelize simulation in x-axis */
    

    if (file.is_open()) {

        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            if (var[0] == '#') { /* ignore comment line*/
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "xlength") file >> xlength;
                if (var == "ylength") file >> ylength;
                if (var == "nu") file >> nu;
                if (var == "t_end") file >> _t_end;
                if (var == "dt") file >> dt;
                if (var == "omg") file >> omg;
                if (var == "eps") file >> eps;
                if (var == "tau") file >> tau;
                if (var == "gamma") file >> gamma;
                if (var == "dt_value") file >> _output_freq;
                if (var == "UI") file >> UI;
                if (var == "VI") file >> VI;
                if (var == "GX") file >> GX;
                if (var == "GY") file >> GY;
                if (var == "PI") file >> PI;
                if (var == "itermax") file >> itermax;
                if (var == "imax") file >> imax;
                if (var == "jmax") file >> jmax;
                if (var == "UIN") file >> UIN;
                if (var == "VIN") file >> VIN;
                if (var == "energy_eq") file >> energy_eq;
	        if (var == "TI") file >> TI;
	        if (var == "beta") file >> beta;
	        if (var == "alpha") file >> alpha;
	        if (var == "num_walls") file >> num_walls;
	        if ((num_walls != 0) and (num_walls != 6)){
	        	 
			for (int i = 3; i < 3 + num_walls; i++){
				if (var == "wall_temp_" + std::to_string(i)) file >> temps[i - 3];
			}
	        }
	        if (var == "geo_file") file >> geo_file;
	        if (var == "iproc") file >> iproc;
	        if (var == "jproc") file >> jproc;
            }
        }
    }
    file.close();

    _datfile_name = file_name;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    _process_rank = process_rank;
    _size = size;

    if (_iproc * _jproc != _size) {
        Communication::finalize();
        if (_process_rank == 0) {
            std::cerr << "please check your iproc(number of processes for x-direction) and jproc(number of processes "
                         "for y-direction). Their product should be the total number of processes for the problem"
                      << std::endl;
        }
        exit(0);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    _iproc = iproc;
    _jproc = jproc;
    
    if ((_iproc == 1) and (_jproc == 1)){
    	std::string question;
	while ((question != "yes") and (question != "no")){
		std::cout << "Specified number of processors on x and y-axis set to 1. Running simulation in sequential mode (only one processor), are you sure you want to continue? Please type 'yes' or 'no' to continue. If using parallelization (more processors) runtime can be significantly reduced:\n";
		std::cin >> question;
	}
	if (question == "no"){
		std::ofstream outFile(file_name, std::ios::app); // Open file in append mode

	        if (!outFile) {
		    std::cerr << "Error opening file " << file_name << " to add iproc and jproc parameters." << std::endl;
		    exit(1);
	        }

		std::cout << "Please input the number of processes for the x-axis:\n";
		std::cin >> iproc;
		
		while (!std::cin.good())
		{
		    std::cin.clear();
		    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		    std::cout << "Please input the number of processes for the x-axis:\n";
		    std::cin >> iproc;
		}
		
		std::cout << "Please input the number of processes for the y-axis:\n";
		std::cin >> jproc;
		while (!std::cin.good())
		{
		    std::cin.clear();
		    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		    std::cout << "Please input the number of processes for the y-axis:\n";
		    std::cin >> jproc;
		}
		std::string lineToAdd1 = "iproc          " + std::to_string(iproc);
		std::string lineToAdd2 = "jproc          " + std::to_string(jproc);

	        outFile << lineToAdd1 << std::endl; // Write the line to the file
	        outFile << lineToAdd2 << std::endl;

	        outFile.close(); // Close the file
	    
	    	std::cout << "Thank you. Successfully added parameters to '.dat' file with " << iproc << " processors in x-axis, and " << jproc << " processors in y-axis." << std::endl;
	    	std::cout << "Proceeding to terminate program. Please run again the scrip using the following command: 'mpirun -np " << iproc*jproc << " ./fluidchen " << file_name << "'" << std::endl;
	    	exit(0);	
	}
	else {
		std::cout << "Proceeding to run simulation with only one processor." << std::endl;
	}
    }

    std::map<int, double> wall_vel;
    std::map<int, double> wall_temp;
    
    if (geo_file.compare("NONE") == 0) {
        wall_vel.insert(std::pair<int, double>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    }
    else {
    	_geom_name = geo_file;
    	if (energy_eq.compare("on") == 0) {
    		_energy_eq = true;
    		for (int i = 3; i < 3 + num_walls; ++i){
    			wall_temp.insert({i, temps[i - 3]}); 
    		}
    	}
    }

    // Set file names for geometry file and output directory
    set_file_names(file_name);

    MPI_Init(&argn, &args);
    //MPI_Init(NULL, NULL);
    int size, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if (_iproc*_jproc != size){
    	std::cout << "Warning: iproc, jproc and number of processes specified at runtime do not match. Terminating program\nPlease make sure to run the code using the following command line: 'mpirun -np " << iproc*jproc << " ./fluidchen " << file_name << "'" << std::endl;
    	exit(1);
    }
    
    std::cout << "iproc: " << _iproc << ". jproc: " << _jproc << ". Total number of threads: " << size << std::endl;
    std::cout << "Total domain should have the following dimentions (number of cells, in x and y: (" << imax << ", " << jmax << ")." << std::endl;
    
    int imin_local;
    int imax_local;
    int jmin_local;
    int jmax_local;
    
    int my_x = my_rank%_iproc;
    int my_y = (my_rank - my_x)/_iproc;
    
    if ((imax % _iproc) == 0){ //All processors take same number of tiles (x-axis).
	imin_local = my_x*(imax/_iproc) + 1;
	imax_local = (my_x + 1)*(imax/_iproc);
	
	if ((jmax % _jproc) == 0){ //All processors take same number of tiles (y-axis).
		jmin_local = my_y*(jmax/_jproc) + 1;
		jmax_local = (my_y + 1)*(jmax/_jproc);
	}
	else{ //Not all processors take same number of tiles (y-axis).
	    	if (my_y < (jmax % _jproc)){ //Processors with an extra tile (y-axis).
	    		jmin_local = my_y*(jmax/_jproc) + (my_y + 1);
			jmax_local = (my_y + 1)*(jmax/_jproc) + (my_y + 1);
	 	}
	 	else{ //Processors with no extra tile (y-axis).
	 		jmin_local = my_y*(jmax/_jproc) + (jmax % _jproc + 1);
			jmax_local = (my_y + 1)*(jmax/_jproc) + (jmax % _jproc);
	 	}
	}
	
	    
    }
    
    else{ //Not all processors take same number of tiles (x-axis).
    	if (my_x < (imax % _iproc)){ //Processors with an extra tile (x-axis).
    		imin_local = my_x*(imax/_iproc) + (my_x + 1);
		imax_local = (my_x + 1)*(imax/_iproc) + (my_x + 1);
 	}
 	else{ //Processors with no extra tile (x-axis).
 		imin_local = my_x*(imax/_iproc) + (imax % _iproc + 1);
		imax_local = (my_x + 1)*(imax/_iproc) + (imax % _iproc);
 	}
 	if ((jmax % _jproc) == 0){ //All processors take same number of tiles (y-axis).
		jmin_local = my_y*(jmax/_jproc) + 1;
		jmax_local = (my_y + 1)*(jmax/_jproc);
	}
	else{ //Not all processors take same number of tiles (y-axis).
	    	if (my_y < (jmax % _jproc)){ //Processors with an extra tile (y-axis).
	    		jmin_local = my_y*(jmax/_jproc) + (my_y + 1);
			jmax_local = (my_y + 1)*(jmax/_jproc) + (my_y + 1);
	 	}
	 	else{ //Processors with no extra tile (y-axis).
	 		jmin_local = my_y*(jmax/_jproc) + (jmax % _jproc + 1);
			jmax_local = (my_y + 1)*(jmax/_jproc) + (jmax % _jproc);
	 	}
	}
    }
    
    std::cout << "I am thread with id: " << my_rank << ". My x and y coordinates are: (" << my_x << ", " << my_y << "). I got assigned tiles from x-position: [" << imin_local << ", " << imax_local << "], and y-position: [" << jmin_local << ", " << jmax_local << "]." << std::endl;
    
    MPI_Finalize();
    
    // Build up the domain
    Domain domain;
    domain.dx = xlength / static_cast<double>(imax);
    domain.dy = ylength / static_cast<double>(jmax);
    domain.domain_size_x = imax;
    domain.domain_size_y = jmax;

    build_domain(domain, imax, jmax);
    
    _grid = Grid(_geom_name, domain);
    _field = Fields(GX, GY, nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, UIN, VIN, PI, TI, alpha, beta);
    Discretization _discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
    _tolerance = eps;	
	
    // Construct boundaries
    if (not _grid.moving_wall_cells().empty()) {
        _boundaries.push_back(
            std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
    }
    if (not _grid.fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells(), wall_temp));
    }
    if (not _grid.inflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<InflowBoundary>(_grid.inflow_cells(), UIN, VIN));
    }
    if (not _grid.outflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<OutflowBoundary>(_grid.outflow_cells()));
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (_process_rank == 0) {
        output_log(_datfile_name, nu, UI, VI, PI, GX, GY, xlength, ylength, dt, imax, jmax, gamma, omg, tau, itermax,
                   eps, TI, alpha, beta, num_walls, my_rank);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (_process_rank == 0) {
    std::cout << "KILL PROGRAM" << std::endl;
}

void Case::set_file_names(std::string file_name) {
    std::string temp_dir;
    bool case_name_flag = true;
    bool prefix_flag = false;

    for (int i = file_name.size() - 1; i > -1; --i) {
        if (file_name[i] == '/') {
            case_name_flag = false;
            prefix_flag = true;
        }
        if (case_name_flag) {
            _case_name.push_back(file_name[i]);
        }
        if (prefix_flag) {
            _prefix.push_back(file_name[i]);
        }
    }

    for (int i = file_name.size() - _case_name.size() - 1; i > -1; --i) {
        temp_dir.push_back(file_name[i]);
    }

    std::reverse(_case_name.begin(), _case_name.end());
    std::reverse(_prefix.begin(), _prefix.end());
    std::reverse(temp_dir.begin(), temp_dir.end());

    _case_name.erase(_case_name.size() - 4);
    _dict_name = temp_dir;
    _dict_name.append(_case_name);
    _dict_name.append("_Output");

    if (_geom_name.compare("NONE") != 0) {
        _geom_name = _prefix + _geom_name;
    }

    // Create output directory
    filesystem::path folder(_dict_name);
    try {
        filesystem::create_directory(folder);
    } catch (const std::exception &e) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}

/**
 * This function is the main simulation loop. In the simulation loop, following steps are required
 * - Calculate and apply boundary conditions for all the boundaries in _boundaries container
 *   using apply() member function of Boundary class
 * - Calculate fluxes (F and G) using calculate_fluxes() member function of Fields class.
 *   Flux consists of diffusion and convection part, which are located in Discretization class
 * - Calculate right-hand-side of PPE using calculate_rs() member function of Fields class
 * - Iterate the pressure poisson equation until the residual becomes smaller than the desired tolerance
 *   or the maximum number of the iterations are performed using solve() member function of PressureSolver
 * - Update boundary conditions after each iteration of the SOR solver
 * - Calculate the velocities u and v using calculate_velocities() member function of Fields class
 * - Calculate the maximal timestep size for the next iteration using calculate_dt() member function of Fields class
 * - Write vtk files using output_vtk() function
 *
 * Please note that some classes such as PressureSolver, Boundary are abstract classes which means they only provide the
 * interface. No member functions should be defined in abstract classes. You need to define functions in inherited
 * classes such as MovingWallBoundary class.
 *
 * For information about the classes and functions, you can check the header files.
 */
void Case::simulate() {

    double t = 0.0;
    double dt = _field.dt();
    int timestep = 0;
    double output_counter = 0.0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     Case::output_vtk(timestep, my_rank);
    std::ofstream output;
    int progress, last_progress;

    if (_process_rank == 0) {
        std::string str = _dict_name + "_run_log_" + std::to_string(my_rank) + ".log";

        output.open(str, std::ios_base::app);
        output << "\n\nIteration Log:\n";
        std::cout << "Simulation is Running!\nPlease Refer to " << _dict_name << "_run_log_" << my_rank
                  << ".log for Simulation log!\n";
        // Simulation Progress
        if (Case::_energy_eq == "on") {
            std::cout << "\nEnergy Equation is On\n";
        }
    }
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
    while (t < _t_end) {
	    
	    dt = _field.calculate_dt(_grid, _energy_eq);
        dt = Communication::reduce_min(dt);
	    
	    if (_energy_eq == true){
	    	_field.calculate_temperatures(_grid, _discretization);
             // communicate temperature
            Communication::communicate(_field.t_matrix(), domain);
	    }
	    
	    _field.calculate_fluxes(_grid, _discretization, _energy_eq);
        // communicate fluxes
        Communication::communicate(_field.f_matrix(), domain);
        Communication::communicate(_field.g_matrix(), domain);
	    for (auto const& boundary : _boundaries){
	    	boundary->apply(_field, _energy_eq);
	    }
	    
	    _field.calculate_rs(_grid);
	    
	    int iter{0};
	    double res = _pressure_solver->solve(_field, _grid, _boundaries); 
	    while ((res > _tolerance) and (iter < _max_iter)) {
	    	res = _pressure_solver->solve(_field, _grid, _boundaries);
            // weighted addition of residuals
            err = Communication::reduce_sum(err);
            fluid_cells = _grid.fluid_cells().size();
            fluid_cells = Communication::reduce_sum(fluid_cells);
            err = std::sqrt(err / fluid_cells);
            // communicate pressures
            Communication::communicate(_field.p_matrix(), domain);
	    	for (auto const& boundary : _boundaries){
	    		boundary->apply(_field, _energy_eq);
	    	}
	    	iter = iter + 1;
	    }
	    
	    if (std::ceil(t) == std::floor(t + dt)){
	    	if (int(std::ceil(t)) % 2 == 0) {
			std::cout << "SOR loop (for calculating pressure at the next time step from the Poisson equation) results:" << '\n';
			std::cout << "t = " << t << ", res = " << res << ", tolerance = " << _tolerance << ", iter = " << iter << ", max iter = " << _max_iter <<  '\n';
		}
	    }
	    
	    _field.calculate_velocities(_grid);
        Communication::communicate(_field.u_matrix(), domain);
        Communication::communicate(_field.v_matrix(), domain);
	    
	    t = t + dt;
    	timestep = timestep + 1;
    	//output_vtk(timestep);
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
        if (_process_rank == 0) {
            output << std::setprecision(4) << std::fixed;
        }
        if (t - output_counter * _output_freq >= 0) {
            Case::output_vtk(timestep, my_rank);
            if (_process_rank == 0) {
                output << "Time: " << t << " Residual: " << err << " PPE Iterations: " << iter_count << std::endl;
                if (iter_count == _max_iter || std::isnan(err)) {
                    std::cout << "The PPE Solver didn't converge for Time = " << t
                              << " Please check the log file and increase max iterations or other parameters for "
                                 "convergence"
                              << "\n";
                }
            }
            output_counter += 1;
        }
        // Printing Simulation Progress
        if (_process_rank == 0) {
            progress = t / t_end * 100;
            if (progress % 10 == 0 && progress != last_progress) {
                std::cout << "Time Step: " << timestep << " Residue: " << err << " PPE Iterations: " << iter_count
                          << std::endl;
                std::cout << "[";
                for (int i = 0; i < progress / 10; i++) {
                    std::cout << "=====";
                }
                if (progress == 100)
                    std::cout << "]";
                else
                    std::cout << ">";
                std::cout << progress << "%\r";
                std::cout.flush();
                last_progress = progress;
            }
        }       
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    }
    if (t_end !=
        (output_counter - 1) * _output_freq) // Recording at t_end if the output frequency is not a multiple of t_end
    {
        Case::output_vtk(timestep, my_rank);
        if (_process_rank == 0) {
            output << "Time Step: " << timestep << " Residue: " << err << " PPE Iterations: " << iter_count
                   << std::endl;
        }
        output_counter += 1;
    }
    if (_process_rank == 0) {
        std::cout << "\nSimulation has ended\n";
    }
    output.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

void Case::output_vtk(int timestep, int my_rank) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = _grid.dx();
    double dy = _grid.dy();
    int i_max = _grid.imax();
    int j_max = _grid.jmax();

    double x = _grid.domain().imin * dx;
    double y = _grid.domain().jmin * dy;

    { y += dy; }
    { x += dx; }

    double z = 0;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().imin * dx;
        { x += dx; }
        for (int row = 0; row < _grid.domain().size_x + 1; row++) {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_grid.domain().size_x + 1, _grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

    // Pressure Array
    vtkSmartPointer<vtkDoubleArray> Pressure = vtkSmartPointer<vtkDoubleArray>::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Velocity Array for cell data
    vtkSmartPointer<vtkDoubleArray> Velocity = vtkSmartPointer<vtkDoubleArray>::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);
    
    // Temperature Array for cell data
    vtkSmartPointer<vtkDoubleArray> Temperature = vtkSmartPointer<vtkDoubleArray>::New();
    if (_energy_eq){
	    
	    Temperature->SetName("temperature");
	    Temperature->SetNumberOfComponents(1);
    }

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0
    double index = 0;

    // Print pressure, velocity and temperature from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double pressure = _field.p(i, j);
            Pressure->InsertNextTuple(&pressure);
            if (_energy_eq){
            	double temperature = _field.t(i, j);
            	Temperature->InsertNextTuple(&temperature);
            }
            vel[0] = (_field.u(i - 1, j) + _field.u(i, j)) * 0.5;
            vel[1] = (_field.v(i, j - 1) + _field.v(i, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
            
            if (_grid.cell(i, j).type() == cell_type::FIXED_WALL){
	        structuredGrid->BlankCell(index);
        }
        index++;
    }
    }

    // Velocity Array for point data
    vtkSmartPointer<vtkDoubleArray> VelocityPoints = vtkSmartPointer<vtkDoubleArray>::New();
    VelocityPoints->SetName("velocity");
    VelocityPoints->SetNumberOfComponents(3);

    // Print Velocity from bottom to top
    for (int j = 0; j < _grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _grid.domain().size_x + 1; i++) {
            vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
            vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
            VelocityPoints->InsertNextTuple(vel);
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetCellData()->AddArray(Velocity);
    structuredGrid->GetPointData()->AddArray(VelocityPoints);
    
    if (_energy_eq){
	    // Add Temperature to Structured Grid
	    structuredGrid->GetCellData()->AddArray(Temperature);
    }

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // // Create Filename
    // std::string outputname =
    //     _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "." + std::to_string(timestep) + ".vtk";
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    // Create Filename
    std::string outputname = _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "_" +
                             std::to_string(_process_rank) + "." + std::to_string(timestep) +
                             ".vtk"; // my_rank is the user's input and _process_rank is the process rank
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {
    	
    domain.imin = 0;
    domain.jmin = 0;
    domain.imax = imax_domain + 2;
    domain.jmax = jmax_domain + 2;
    domain.size_x = imax_domain;
    domain.size_y = jmax_domain;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void Case::output_log(std::string dat_file_name, double nu, double UI, double VI, double PI, double GX, double GY,
                      double xlength, double ylength, double dt, double imax, double jmax, double gamma, double omg,
                      double tau, double itermax, double eps, double TI, double alpha, double beta, double num_walls,
                      int my_rank) {

    std::string str = _dict_name + "_run_log_" + std::to_string(my_rank) + ".log";
    std::stringstream stream;
    std::ofstream output(str);

    output << "Log File for : " << dat_file_name << "\n";
    output << "Simulation Parameters:\n";
    output << "xlength : " << xlength << "\n";
    output << "ylength : " << ylength << "\n";
    output << "nu : " << nu << "\n";
    output << "t_end : " << _t_end << "\n";
    output << "dt : " << dt << "\n";
    output << "omg : " << omg << "\n";
    output << "eps : " << eps << "\n";
    output << "tau : " << tau << "\n";
    output << "gamma : " << gamma << "\n";
    output << "dt_value : " << _output_freq << "\n";
    output << "UI : " << UI << "\n";
    output << "VI : " << VI << "\n";
    output << "GX : " << GX << "\n";
    output << "GY : " << GY << "\n";
    output << "PI : " << PI << "\n";
    output << "itermax : " << itermax << "\n";
    output << "Energy Equation : " << _energy_eq << "\n";
    output << "Number of processes in x-direction : " << _iproc << "\n";
    output << "Number of processes in y-direction : " << _jproc << "\n";
    output << "Total number of process : " << _size << "\n";
    if (_energy_eq == "on") {
        output << "Temp Initial : " << TI << "\n";
        output << "alpha : " << alpha << "\n";
        output << "beta : " << beta << "\n";
        output << "No of Temperature walls: " << num_walls << "\n";
        
    }

    output.close();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////