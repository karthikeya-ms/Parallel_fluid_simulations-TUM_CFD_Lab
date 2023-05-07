#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {
    int i_idx{0};
	int j_idx{0};

	for (const auto cell : _cells){
	
		i_idx = cell.i();
		j_idx = cell.j();
        bool tw = cell.is_border(border_position::TOP);
        bool rw = cell.is_border(border_position::RIGHT); 
        bool bw = cell.is_border(border_position::BOTTOM);
        bool lw = cell.is_border(border_position::LEFT);
        bool con = true;

        switch(con){
            case tw:
                field.p(i_idx, j_idx) = field.p(i_idx, j_idx + 1);
			    field.u(i_idx, j_idx) = -field.u(i_idx, j_idx + 1);
			    field.v(i_idx, j_idx) = 0; 
			    field.g(i_idx, j_idx) = field.v(i_idx, j_idx);
                break;
            case rw:
                field.p(i_idx, j_idx) = field.p(i_idx + 1, j_idx);
			    field.u(i_idx, j_idx) = 0;
			    field.v(i_idx, j_idx) = -field.v(i_idx + 1, j_idx); 
			    field.f(i_idx, j_idx) = field.u(i_idx, j_idx);
                break;
            case bw:
                field.p(i_idx, j_idx) = field.p(i_idx, j_idx - 1);
			    field.u(i_idx, j_idx) = -field.u(i_idx, j_idx - 1);
			    field.v(i_idx, j_idx - 1) = 0; 
			    field.g(i_idx, j_idx - 1) = field.v(i_idx, j_idx - 1);
                break;
            case lw:
                field.p(i_idx, j_idx) = field.p(i_idx - 1, j_idx);
			    field.u(i_idx - 1, j_idx) = 0;
			    field.v(i_idx, j_idx) = -field.v(i_idx - 1, j_idx); 
			    field.f(i_idx - 1, j_idx) = field.u(i_idx - 1, j_idx);
                break;

        }
		
	}

}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {

	int i_idx{0};
	int j_idx{0};
	double vel = _wall_velocity.at(LidDrivenCavity::moving_wall_id);

	for (const auto cell : _cells){
	
		i_idx = cell.i();
		j_idx = cell.j();
        bool tw = cell.is_border(border_position::TOP);
        bool rw = cell.is_border(border_position::RIGHT); 
        bool bw = cell.is_border(border_position::BOTTOM);
        bool lw = cell.is_border(border_position::LEFT);
        bool con = true;
		
	    switch(con){
            case tw:
                field.p(i_idx, j_idx) = field.p(i_idx, j_idx + 1);
			    field.u(i_idx, j_idx) = 2*vel - field.u(i_idx, j_idx + 1);
			    field.v(i_idx, j_idx) = 0; 
			    field.g(i_idx, j_idx) = field.v(i_idx, j_idx);
                break;
            case rw:
                field.p(i_idx, j_idx) = field.p(i_idx + 1, j_idx);
			    field.u(i_idx, j_idx) = vel;
			    field.v(i_idx, j_idx) = -field.v(i_idx + 1, j_idx); 
			    field.f(i_idx, j_idx) = field.u(i_idx, j_idx);
                break;
            case bw:
                field.p(i_idx, j_idx) = field.p(i_idx, j_idx - 1);
			    field.u(i_idx, j_idx) = 2*vel - field.u(i_idx, j_idx - 1);
			    field.v(i_idx, j_idx - 1) = 0; 
			    field.g(i_idx, j_idx - 1) = field.v(i_idx, j_idx - 1);
                break;
            case lw:
                field.p(i_idx, j_idx) = field.p(i_idx - 1, j_idx);
			    field.u(i_idx - 1, j_idx) = vel;
			    field.v(i_idx, j_idx) = -field.v(i_idx - 1, j_idx); 
			    field.f(i_idx - 1, j_idx) = field.u(i_idx - 1, j_idx);
                break;

        }
	}

}
