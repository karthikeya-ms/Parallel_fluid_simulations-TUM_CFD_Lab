#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

//Added: Function implemented for first task.
void FixedWallBoundary::apply(Fields &field, bool energy_eq) {

	int i_idx{0};
	int j_idx{0};

	for (const auto currentCell : _cells){
	
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		
		//ASSUMING GIVEN LABELING OF WALLS: 3) ADIABATIC
		//				    4) HOT
		//				    5) COLD
		//FOR DIRICHLET AND ADIABATIC WALLS WITH TWO BORDERS, WE CONTROL TEMPERATURE ON BOTH WALLS AT THE SAME TIME WITH AVERAGES OF BOTH NEIGHBOORS' T.
		
		int num_borders = 0;
		for (auto i = border_position::TOP; i != border_position::LAST; i = static_cast<border_position>((size_t)i + 1)){
			if (currentCell->is_border(i)){
				++num_borders;
			}
		}
		
		if (num_borders == 1){
			if (currentCell->is_border(border_position::TOP) == true){
				field.p(i_idx, j_idx) = field.p(i_idx, j_idx + 1);
				field.u(i_idx, j_idx) = -field.u(i_idx, j_idx + 1);
				field.v(i_idx, j_idx) = 0; 
				field.g(i_idx, j_idx) = field.v(i_idx, j_idx);
				
				if (energy_eq == true){
					if (currentCell->wall_id() == 3){ //Adiabatic walls: Neumann boundary condition for T, with heat flux q_N = 0
						field.t(i_idx, j_idx) = field.t(i_idx, j_idx + 1);
					}
					else if ((currentCell->wall_id() == 4) || (currentCell->wall_id() == 5)){ //Hot/cold walls: Dirichlet boundary condition for T
						field.t(i_idx, j_idx) = 2*_wall_temperature.at(currentCell->wall_id()) - field.t(i_idx, j_idx + 1);
					}
				}
			}
			else if (currentCell->is_border(border_position::BOTTOM) == true) {
				field.p(i_idx, j_idx) = field.p(i_idx, j_idx - 1);
				field.u(i_idx, j_idx) = -field.u(i_idx, j_idx - 1);
				field.v(i_idx, j_idx - 1) = 0; 
				field.g(i_idx, j_idx - 1) = field.v(i_idx, j_idx - 1);
				
				if (energy_eq == true){
					if (currentCell->wall_id() == 3){ //Adiabatic walls: Neumann boundary condition for T, with heat flux q_N = 0
						field.t(i_idx, j_idx) = field.t(i_idx, j_idx - 1);
					}
					else if ((currentCell->wall_id() == 4) || (currentCell->wall_id() == 5)){ //Hot/cold walls: Dirichlet boundary condition for T
						field.t(i_idx, j_idx) = 2*_wall_temperature.at(currentCell->wall_id()) - field.t(i_idx, j_idx - 1);
					}
				}
			}
			else if (currentCell->is_border(border_position::RIGHT) == true) {
				field.p(i_idx, j_idx) = field.p(i_idx + 1, j_idx);
				field.u(i_idx, j_idx) = 0;
				field.v(i_idx, j_idx) = -field.v(i_idx + 1, j_idx); 
				field.f(i_idx, j_idx) = field.u(i_idx, j_idx);
				
				if (energy_eq == true){
					if (currentCell->wall_id() == 3){ //Adiabatic walls: Neumann boundary condition for T, with heat flux q_N = 0
						field.t(i_idx, j_idx) = field.t(i_idx + 1, j_idx);
					}
					else if ((currentCell->wall_id() == 4) || (currentCell->wall_id() == 5)){ //Hot/cold walls: Dirichlet boundary condition for T
						field.t(i_idx, j_idx) = 2*_wall_temperature.at(currentCell->wall_id()) - field.t(i_idx + 1, j_idx);
					}
				}
			}
			else if (currentCell->is_border(border_position::LEFT) == true) {
				field.p(i_idx, j_idx) = field.p(i_idx - 1, j_idx);
				field.u(i_idx - 1, j_idx) = 0;
				field.v(i_idx, j_idx) = -field.v(i_idx - 1, j_idx); 
				field.f(i_idx - 1, j_idx) = field.u(i_idx - 1, j_idx);
				
				if (energy_eq == true){
					if (currentCell->wall_id() == 3){ //Adiabatic walls: Neumann boundary condition for T, with heat flux q_N = 0
						field.t(i_idx, j_idx) = field.t(i_idx - 1, j_idx);
					}
					else if ((currentCell->wall_id() == 4) || (currentCell->wall_id() == 5)){ //Hot/cold walls: Dirichlet boundary condition for T
						field.t(i_idx, j_idx) = 2*_wall_temperature.at(currentCell->wall_id()) - field.t(i_idx - 1, j_idx);
					}
				}
			}
		}
		
		else if (num_borders == 2){
			if ((currentCell->is_border(border_position::TOP) == true) && (currentCell->is_border(border_position::RIGHT) == true)){
				field.p(i_idx, j_idx) = (field.p(i_idx + 1, j_idx) + field.p(i_idx, j_idx + 1))/2.0;
				field.u(i_idx, j_idx) = 0;
				field.u(i_idx - 1, j_idx) = -field.u(i_idx - 1, j_idx + 1);
				field.v(i_idx, j_idx) = 0;
				field.v(i_idx, j_idx - 1) = -field.v(i_idx + 1, j_idx - 1);
				field.f(i_idx, j_idx) = field.u(i_idx, j_idx);
				field.g(i_idx, j_idx) = field.v(i_idx, j_idx);
				
				if (energy_eq == true){
					if (currentCell->wall_id() == 3){ //Adiabatic walls: Neumann boundary condition for T, with heat flux q_N = 0
						field.t(i_idx, j_idx) = (field.t(i_idx, j_idx + 1) + field.t(i_idx + 1, j_idx))/2.0;
					}
					else if ((currentCell->wall_id() == 4) || (currentCell->wall_id() == 5)){ //Hot/cold walls: Dirichlet boundary condition for T
						field.t(i_idx, j_idx) = 2*_wall_temperature.at(currentCell->wall_id()) - (field.t(i_idx, j_idx + 1) + field.t(i_idx + 1, j_idx))/2.0;
					}
				}
			}
			else if ((currentCell->is_border(border_position::TOP) == true) && (currentCell->is_border(border_position::LEFT) == true)) {
				field.p(i_idx, j_idx) = (field.p(i_idx, j_idx + 1) + field.p(i_idx - 1, j_idx))/2.0;
				field.u(i_idx - 1, j_idx) = 0;
				field.u(i_idx, j_idx) = -field.u(i_idx, j_idx + 1);
				field.v(i_idx, j_idx - 1) = -field.v(i_idx - 1, j_idx - 1);
				field.v(i_idx, j_idx) = 0;
				field.f(i_idx - 1, j_idx) = field.u(i_idx - 1, j_idx);
				field.g(i_idx, j_idx) = field.v(i_idx, j_idx);
				
				if (energy_eq == true){
					if (currentCell->wall_id() == 3){ //Adiabatic walls: Neumann boundary condition for T, with heat flux q_N = 0
						field.t(i_idx, j_idx) = (field.t(i_idx, j_idx + 1) + field.t(i_idx - 1, j_idx))/2.0;
					}
					else if ((currentCell->wall_id() == 4) || (currentCell->wall_id() == 5)){ //Hot/cold walls: Dirichlet boundary condition for T
						field.t(i_idx, j_idx) = 2*_wall_temperature.at(currentCell->wall_id()) - (field.t(i_idx, j_idx + 1) + field.t(i_idx - 1, j_idx))/2.0;
					}
				}
			}
			else if ((currentCell->is_border(border_position::BOTTOM) == true) && (currentCell->is_border(border_position::RIGHT) == true)) {
				field.p(i_idx, j_idx) = (field.p(i_idx + 1, j_idx) + field.p(i_idx, j_idx - 1))/2.0;
				field.u(i_idx, j_idx) = 0;
				field.u(i_idx - 1, j_idx) = field.u(i_idx - 1, j_idx - 1);
				field.v(i_idx, j_idx) = -field.v(i_idx + 1, j_idx);
				field.v(i_idx, j_idx - 1) = 0; 
				field.f(i_idx, j_idx) = field.u(i_idx, j_idx);
				field.g(i_idx, j_idx - 1) = field.v(i_idx, j_idx - 1);
				
				if (energy_eq == true){
					if (currentCell->wall_id() == 3){ //Adiabatic walls: Neumann boundary condition for T, with heat flux q_N = 0
						field.t(i_idx, j_idx) = (field.t(i_idx, j_idx - 1) + field.t(i_idx + 1, j_idx))/2.0;
					}
					else if ((currentCell->wall_id() == 4) || (currentCell->wall_id() == 5)){ //Hot/cold walls: Dirichlet boundary condition for T
						field.t(i_idx, j_idx) = 2*_wall_temperature.at(currentCell->wall_id()) - (field.t(i_idx, j_idx - 1) + field.t(i_idx + 1, j_idx))/2.0;
					}
				}
			}
			else if ((currentCell->is_border(border_position::BOTTOM) == true) && (currentCell->is_border(border_position::LEFT) == true)) {
				field.p(i_idx, j_idx) = (field.p(i_idx - 1, j_idx) + field.p(i_idx, j_idx - 1))/2.0;
				field.u(i_idx, j_idx) = 0;
				field.u(i_idx - 1, j_idx) = -field.u(i_idx - 1, j_idx - 1);
				field.v(i_idx, j_idx) = -field.v(i_idx + 1, j_idx);
				field.v(i_idx, j_idx - 1) = 0; 
				field.f(i_idx, j_idx) = field.u(i_idx, j_idx);
				field.g(i_idx, j_idx - 1) = field.v(i_idx, j_idx - 1);
				
				if (energy_eq == true){
					if (currentCell->wall_id() == 3){ //Adiabatic walls: Neumann boundary condition for T, with heat flux q_N = 0
						field.t(i_idx, j_idx) = (field.t(i_idx, j_idx - 1) + field.t(i_idx - 1, j_idx))/2.0;
					}
					else if ((currentCell->wall_id() == 4) || (currentCell->wall_id() == 5)){ //Hot/cold walls: Dirichlet boundary condition for T
						field.t(i_idx, j_idx) = 2*_wall_temperature.at(currentCell->wall_id()) - (field.t(i_idx, j_idx - 1) + field.t(i_idx - 1, j_idx))/2.0;
					}
				}
			}
		}
		
		else if (num_borders > 2){
			std::cout << "Incorrect '.pgm' geometry file. Border cells cannot have more than two neighboors: cell at position (" << i_idx << ", " << j_idx << ") has more than two neighboors. Please correct '.pgm' file and compile again.\n";
		}
	}
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

//Added: Function implemented for first task.
void MovingWallBoundary::apply(Fields &field, bool energy_eq) {

	int i_idx{0};
	int j_idx{0};
	double vel = _wall_velocity.at(LidDrivenCavity::moving_wall_id);

	for (const auto currentCell : _cells){
	
		i_idx = currentCell->i();
		j_idx = currentCell->j();
		
		if (currentCell->is_border(border_position::TOP) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx, j_idx + 1);
			field.u(i_idx, j_idx) = 2*vel - field.u(i_idx, j_idx + 1);
			field.v(i_idx, j_idx) = 0; 
			field.g(i_idx, j_idx) = field.v(i_idx, j_idx);
		}
		else if (currentCell->is_border(border_position::LEFT) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx - 1, j_idx);
			field.u(i_idx - 1, j_idx) = vel;
			field.v(i_idx, j_idx) = -field.v(i_idx - 1, j_idx); 
			field.f(i_idx - 1, j_idx) = field.u(i_idx - 1, j_idx);
		}
		else if (currentCell->is_border(border_position::BOTTOM) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx, j_idx - 1);
			field.u(i_idx, j_idx) = 2*vel - field.u(i_idx, j_idx - 1);
			field.v(i_idx, j_idx - 1) = 0; 
			field.g(i_idx, j_idx - 1) = field.v(i_idx, j_idx - 1);
		}
		else if (currentCell->is_border(border_position::RIGHT) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx + 1, j_idx);
			field.u(i_idx, j_idx) = vel;
			field.v(i_idx, j_idx) = -field.v(i_idx + 1, j_idx); 
			field.f(i_idx, j_idx) = field.u(i_idx, j_idx);
		}
	}
}

InflowBoundary::InflowBoundary(std::vector<Cell *> cells) : _cells(cells) {}
InflowBoundary::InflowBoundary(std::vector<Cell *> cells, double UIN, double VIN) : _cells(cells), _UIN(UIN), _VIN(VIN) {}

void InflowBoundary::apply(Fields &field, bool energy_eq) {
	
	int i_idx{0};
	int j_idx{0};
	for (const auto currentCell : _cells){
	
		i_idx = currentCell->i();
		j_idx = currentCell->j();
	
		// field.p(i_idx, j_idx) = field.p(i_idx + 1, j_idx);
		// field.u(i_idx, j_idx) = _UIN;
       	// field.v(i_idx, j_idx) = 2 * _VIN - field.v(i_idx + 1, j_idx); 
		// field.f(i_idx, j_idx) = field.u(i_idx, j_idx);

		if (currentCell->is_border(border_position::TOP) == true) {
			field.u(i_idx, j_idx) = _UIN; // U component velocity at the top border
			field.v(i_idx, j_idx) = _VIN; // V component velocity at the top border
			field.p(i_idx, j_idx) = field.p(i_idx, j_idx + 1); // Zero Neumann condition for pressure
			field.g(i_idx, j_idx) = field.v(i_idx, j_idx); // Dirichlet condition for V component flux
		}
		else if (currentCell->is_border(border_position::LEFT) == true) {
			field.u(i_idx - 1, j_idx) = _UIN; // U component velocity at the left border
			field.v(i_idx, j_idx) = _VIN; // V component velocity at the left border
			field.p(i_idx, j_idx) = field.p(i_idx - 1, j_idx); // Zero Neumann condition for pressure
			field.f(i_idx - 1, j_idx) = field.u(i_idx - 1, j_idx); // Dirichlet condition for U component flux
		}
		else if (currentCell->is_border(border_position::BOTTOM) == true) {
			field.u(i_idx, j_idx) = _UIN; // U component velocity at the bottom border
			field.v(i_idx, j_idx - 1) = _VIN; // V component velocity at the bottom border
			field.p(i_idx, j_idx) = field.p(i_idx, j_idx - 1); // Zero Neumann condition for pressure
			field.g(i_idx, j_idx - 1) = field.v(i_idx, j_idx - 1); // Dirichlet condition for V component flux
		}
		else if (currentCell->is_border(border_position::RIGHT) == true) {
			field.u(i_idx, j_idx) = _UIN; // U component velocity at the right border
			field.v(i_idx + 1, j_idx) = _VIN; // V component velocity at the right border
			field.p(i_idx, j_idx) = field.p(i_idx + 1, j_idx); // Zero Neumann condition for pressure
			field.f(i_idx, j_idx) = field.u(i_idx, j_idx); // Dirichlet condition for U component flux
		}
	}
}

OutflowBoundary::OutflowBoundary(std::vector<Cell *> cells) : _cells(cells) {}

void OutflowBoundary::apply(Fields &field, bool energy_eq) {

	int i_idx{0};
	int j_idx{0};
	for (const auto currentCell : _cells){
	
		// field.p(i_idx, j_idx) = field.p(i_idx - 1, j_idx);
		// field.u(i_idx, j_idx) = field.u(i_idx - 1, j_idx);
		// field.v(i_idx, j_idx) = field.v(i_idx - 1, j_idx);
		// field.v(i_idx, j_idx - 1) = field.v(i_idx - 1, j_idx - 1); 
		// field.f(i_idx - 1, j_idx) = field.u(i_idx, j_idx);


		if (currentCell->is_border(border_position::TOP) == true) {
			field.p(i_idx, j_idx) = 0;
			field.u(i_idx, j_idx) = field.u(i_idx, j_idx + 1);
			field.v(i_idx, j_idx) = field.v(i_idx, j_idx + 1);
			field.g(i_idx, j_idx) = field.v(i_idx, j_idx);
		}
		else if (currentCell->is_border(border_position::LEFT) == true) {
			field.p(i_idx, j_idx) = 0;
			field.u(i_idx - 1, j_idx) = field.u(i_idx, j_idx);
			field.v(i_idx, j_idx) = field.v(i_idx - 1, j_idx);
			field.f(i_idx - 1, j_idx) = field.u(i_idx - 1, j_idx);
		}
		else if (currentCell->is_border(border_position::BOTTOM) == true) {
			field.p(i_idx, j_idx) = 0;
			field.u(i_idx, j_idx) = field.u(i_idx, j_idx - 1);
			field.v(i_idx, j_idx) = field.v(i_idx, j_idx - 1);
			field.g(i_idx, j_idx - 1) = field.v(i_idx, j_idx - 1);
		}
		else if (currentCell->is_border(border_position::RIGHT) == true) {
			field.p(i_idx, j_idx) = 0;
			field.u(i_idx, j_idx) = field.u(i_idx + 1, j_idx);
            field.v(i_idx + 1, j_idx) = field.v(i_idx, j_idx);
			field.f(i_idx, j_idx) = field.u(i_idx, j_idx);
		}
	}
}
