#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

//Added: Function implemented for first task.
void FixedWallBoundary::apply(Fields &field) {

	int i_idx{0};
	int j_idx{0};

	for (const auto cell : _cells){
	
		i_idx = cell->i();
		j_idx = cell->j();
		
		if (cell->is_border(border_position::TOP) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx, j_idx + 1);
			field.u(i_idx, j_idx) = -field.u(i_idx, j_idx + 1);
			field.v(i_idx, j_idx) = 0; 
			field.g(i_idx, j_idx) = field.v(i_idx, j_idx);
		}
		else if (cell->is_border(border_position::LEFT) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx - 1, j_idx);
			field.u(i_idx - 1, j_idx) = 0;
			field.v(i_idx, j_idx) = -field.v(i_idx - 1, j_idx); 
			field.f(i_idx - 1, j_idx) = field.u(i_idx - 1, j_idx);
		}
		else if (cell->is_border(border_position::BOTTOM) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx, j_idx - 1);
			field.u(i_idx, j_idx) = -field.u(i_idx, j_idx - 1);
			field.v(i_idx, j_idx - 1) = 0; 
			field.g(i_idx, j_idx - 1) = field.v(i_idx, j_idx - 1);
		}
		else if (cell->is_border(border_position::RIGHT) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx + 1, j_idx);
			field.u(i_idx, j_idx) = 0;
			field.v(i_idx, j_idx) = -field.v(i_idx + 1, j_idx); 
			field.f(i_idx, j_idx) = field.u(i_idx, j_idx);
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
void MovingWallBoundary::apply(Fields &field) {

	int i_idx{0};
	int j_idx{0};
	double vel = _wall_velocity.at(LidDrivenCavity::moving_wall_id);

	for (const auto cell : _cells){
	
		i_idx = cell->i();
		j_idx = cell->j();
		
		if (cell->is_border(border_position::TOP) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx, j_idx + 1);
			field.u(i_idx, j_idx) = 2*vel - field.u(i_idx, j_idx + 1);
			field.v(i_idx, j_idx) = 0; 
			field.g(i_idx, j_idx) = field.v(i_idx, j_idx);
		}
		else if (cell->is_border(border_position::LEFT) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx - 1, j_idx);
			field.u(i_idx - 1, j_idx) = vel;
			field.v(i_idx, j_idx) = -field.v(i_idx - 1, j_idx); 
			field.f(i_idx - 1, j_idx) = field.u(i_idx - 1, j_idx);
		}
		else if (cell->is_border(border_position::BOTTOM) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx, j_idx - 1);
			field.u(i_idx, j_idx) = 2*vel - field.u(i_idx, j_idx - 1);
			field.v(i_idx, j_idx - 1) = 0; 
			field.g(i_idx, j_idx - 1) = field.v(i_idx, j_idx - 1);
		}
		else if (cell->is_border(border_position::RIGHT) == true) {
			field.p(i_idx, j_idx) = field.p(i_idx + 1, j_idx);
			field.u(i_idx, j_idx) = vel;
			field.v(i_idx, j_idx) = -field.v(i_idx + 1, j_idx); 
			field.f(i_idx, j_idx) = field.u(i_idx, j_idx);
		}
	}
}
