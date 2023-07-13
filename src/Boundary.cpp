#include "Boundary.hpp"
#include <cmath>
#include <iostream>


FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, double wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {
    }

void FixedWallBoundary::apply(Fields &field) {
    for(const auto& currentCell: _cells){
        
        if(currentCell->is_border(border_position::TOP)){ //Fluid on top cells
            // V(i,j) = 0.0
            field.setv(currentCell->i(),currentCell->j(),0.0);     
            // U(i,j) = -U(i,j+1)
            field.setu(currentCell->i(),currentCell->j(),-(field.u(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j())));
            // P(i,j) = P(i,j+1)
            field.setp(currentCell->i(),currentCell->j(),field.p(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()));
            // G(i,j) = V(i,j)
            field.setg(currentCell->i(),currentCell->j(),field.v(currentCell->i(),currentCell->j()));
            if(field.Energy() == "on"){
                if (currentCell->type() == cell_type::ADIABATIC_WALL){
                    // T(i,j) = T(i,j+1)
                    field.setT(currentCell->i(),currentCell->j(),field.T(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()));            
                }
                if (currentCell->type() == cell_type::HOT_WALL || currentCell->type() == cell_type::COLD_WALL){
                    // T(i,j) = 2*T_wall - T(i,j+1)
                    field.setT(currentCell->i(),currentCell->j(),2*_wall_temperature - field.T(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()));
                }  
            }
        }
        
        if(currentCell->is_border(border_position::BOTTOM)){ //Fluid on bottom cells
            // V(i,j) = 0.0
            field.setv(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),0.0);     
            // U(i,j) = -U(i,j-1)
            field.setu(currentCell->i(),currentCell->j(),-(field.u(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j())));
            // P(i,j) = P(i,j-1)
            field.setp(currentCell->i(),currentCell->j(),field.p(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
            // G(i,j-1) = V(i,j-1)
            field.setg(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),field.v(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
            if(field.Energy() == "on"){
                if (currentCell->type() == cell_type::ADIABATIC_WALL){
                    // T(i,j) = T(i,j-1)
                    field.setT(currentCell->i(),currentCell->j(),field.T(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));            
                }
                if (currentCell->type() == cell_type::HOT_WALL || currentCell->type() == cell_type::COLD_WALL){
                    // T(i,j) = 2*T_wall - T(i,j-1)
                    field.setT(currentCell->i(),currentCell->j(),2*_wall_temperature - field.T(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
                }  
            }        
        }

        if(currentCell->is_border(border_position::RIGHT)){   //Fluid on right cells,
            // U(i,j) = 0.0
            field.setu(currentCell->i(),currentCell->j(),0.0);
            // V(i,j) = -V(i+1,j)
            field.setv(currentCell->i(),currentCell->j(),-(field.v(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j())));
            // P(i,j) = P(i+1,j)
            field.setp(currentCell->i(),currentCell->j(),field.p(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j()));
            // F(i,j) = U(i,j)
            field.setf(currentCell->i(),currentCell->j(),field.u(currentCell->i(),currentCell->j()));
            if(field.Energy() == "on"){
                if (currentCell->type() == cell_type::ADIABATIC_WALL){
                    // T(i,j) = T(i+1,j)
                    field.setT(currentCell->i(),currentCell->j(),field.T(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j()));            
                }
                if (currentCell->type() == cell_type::HOT_WALL || currentCell->type() == cell_type::COLD_WALL){
                    // T(i,j) = 2*T_wall - T(i+1,j)
                    field.setT(currentCell->i(),currentCell->j(),2*_wall_temperature - field.T(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j()));
                }  
            }   
        }       
        if(currentCell->is_border(border_position::LEFT)){    //Fluid on left cells
            // U(i,j) = 0.0
            field.setu(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),0.0);
            // V(i,j) = -V(i-1,j)
            field.setv(currentCell->i(),currentCell->j(),-(field.v(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j())));
            // P(i,j) = P(i-1,j)
            field.setp(currentCell->i(),currentCell->j(),field.p(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));
            // F(i-1,j) = U(i-1,j)
            field.setf(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),field.u(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));
            if(field.Energy() == "on"){
                if (currentCell->type() == cell_type::ADIABATIC_WALL){
                    // T(i,j) = T(i-1,j)
                    field.setT(currentCell->i(),currentCell->j(),field.T(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));            
                }
                if (currentCell->type() == cell_type::HOT_WALL || currentCell->type() == cell_type::COLD_WALL){
                    // T(i,j) = 2*T_wall - T(i-1,j)
                    field.setT(currentCell->i(),currentCell->j(),2*_wall_temperature - field.T(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));
                }  
            }  
        }

        if(currentCell->is_border(border_position::TOP) && currentCell->is_border(border_position::RIGHT)){ //Fluid on top and on right cells
            // V(i,j) = 0
            field.setv(currentCell->i(),currentCell->j(),0.0);    
            // U(i,j) = 0
            field.setu(currentCell->i(),currentCell->j(),0.0);     
            // U(i-1,j) = -U(i-1,j+1)
            field.setu(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),-field.u(currentCell->neighbour(border_position::LEFT)->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::LEFT)->neighbour(border_position::TOP)->j()));  
            // V(i,j-1) = -V(i+1,j-1)
            field.setv(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),-field.v(currentCell->neighbour(border_position::BOTTOM)->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::BOTTOM)->neighbour(border_position::RIGHT)->j()));  
            //P_(i,j) = (P_(i,j+1)+P_(i+1,j))/2
            field.setp(currentCell->i(),currentCell->j(),(field.p(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()) + field.p(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j()))*0.5);  
            // F(i,j) = U(i,j)
            field.setf(currentCell->i(),currentCell->j(),field.u(currentCell->i(),currentCell->j()));
            // G(i,j) = V(i,j)
            field.setg(currentCell->i(),currentCell->j(),field.v(currentCell->i(),currentCell->j()));
            if(field.Energy() == "on"){
                if (currentCell->type() == cell_type::ADIABATIC_WALL){
                    // T(i,j) = (T(i,j+1) + T(i+1,j))/2
                    field.setT(currentCell->i(),currentCell->j(),0.5*(field.T(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()) + field.T(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j())));            
                }
                if (currentCell->type() == cell_type::HOT_WALL || currentCell->type() == cell_type::COLD_WALL){
                    // T(i,j) = 2*T_wall - (T(i,j+1) + T(i+1,j))/2
                    field.setT(currentCell->i(),currentCell->j(),2*_wall_temperature - 0.5* (field.T(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()) + field.T(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j())));
                }  
            }  
        }

        if(currentCell->is_border(border_position::TOP) && currentCell->is_border(border_position::LEFT)){ //Fluid on top and on left cells
            // V(i,j) = 0.0
            field.setv(currentCell->i(),currentCell->j(),0.0);    
            // U(i-1,j) = 0.0
            field.setu(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),0.0);     
            // U(i,j) = -U(i,j+1)
            field.setu(currentCell->i(),currentCell->j(),-field.u(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()));  
            // V(i,j-1) = -V(i-1,j-1)
            field.setv(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),-field.v(currentCell->neighbour(border_position::BOTTOM)->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::BOTTOM)->neighbour(border_position::LEFT)->j()));  
            // P(i,j) = (P(i,j+1)+P(i-1,j))/2
            field.setp(currentCell->i(),currentCell->j(),(field.p(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()) + field.p(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()))*0.5);  
            // F(i-1,j) = U(i-1,j)
            field.setf(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),field.u(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));
            // G(i,j) = V(i,j)
            field.setg(currentCell->i(),currentCell->j(),field.v(currentCell->i(),currentCell->j()));
            if(field.Energy() == "on"){
                if (currentCell->type() == cell_type::ADIABATIC_WALL){
                    // T(i,j) = (T(i,j+1)+T(i-1,j))/2
                    field.setT(currentCell->i(),currentCell->j(),0.5*(field.T(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()) + field.T(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j())));            
                }
                if (currentCell->type() == cell_type::HOT_WALL || currentCell->type() == cell_type::COLD_WALL){
                    // T(i,j) = 2*T_wall - (T(i,j+1)+T(i-1,j))/2
                    field.setT(currentCell->i(),currentCell->j(),2*_wall_temperature - 0.5* (field.T(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()) + field.T(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j())));
                } 
            } 
        }

        if(currentCell->is_border(border_position::BOTTOM) && currentCell->is_border(border_position::RIGHT)){ //Fluid on bottom and on right cells
            // V(i,j-1) = 0.0
            field.setv(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),0.0);     
            // U(i,j) = 0.0
            field.setu(currentCell->i(),currentCell->j(),0.0);     
            // U(i-1,j) = -U(i-1,j-1)
            field.setu(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),-field.u(currentCell->neighbour(border_position::LEFT)->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::LEFT)->neighbour(border_position::BOTTOM)->j()));  
            //V(i,j) = -V(i+1,j)
            field.setv(currentCell->i(),currentCell->j(),-field.v(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j()));  
            // P(i,j) = (P(i,j-1)+P(i+1,j))/2
            field.setp(currentCell->i(),currentCell->j(),(field.p(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()) + field.p(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j()))*0.5);  
            // F(i,j) = U(i,j)
            field.setf(currentCell->i(),currentCell->j(),field.u(currentCell->i(),currentCell->j()));
            // G(i,j-1) = V(i,j-1)
            field.setg(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),field.v(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
            if(field.Energy() == "on"){
                if (currentCell->type() == cell_type::ADIABATIC_WALL){
                    // T(i,j) = (T(i,j-1)+T(i+1,j))/2
                    field.setT(currentCell->i(),currentCell->j(),0.5*(field.T(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()) + field.T(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j())));            
                }
                if (currentCell->type() == cell_type::HOT_WALL || currentCell->type() == cell_type::COLD_WALL){
                    // T(i,j) = 2*T_wall - (T(i,j-1)+T(i+1,j))/2
                    field.setT(currentCell->i(),currentCell->j(),2*_wall_temperature - 0.5* (field.T(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()) + field.T(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j())));
                } 
            } 
        }

        if(currentCell->is_border(border_position::BOTTOM) && currentCell->is_border(border_position::LEFT)){ //Fluid on bottom and on left cells
            // V(i,j-1) = 0.0
            field.setv(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),0.0);     
            // U(i-1,j-1) = 0.0
            field.setu(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),0.0);     
            // U(i,j) = -U(i,j-1)
            field.setu(currentCell->i(),currentCell->j(),-field.u(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));  
            // V(i,j) = -V(i-1,j)
            field.setv(currentCell->i(),currentCell->j(),-field.v(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));  
            // P(i,j) = (P(i,j-1)+P(i-1,j))/2
            field.setp(currentCell->i(),currentCell->j(),(field.p(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()) + field.p(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()))*0.5);  
            // F(i-1,j) = U(i-1,j)
            field.setf(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),field.u(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));
            // G(i,j-1) = V(i,j-1)
            field.setg(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),field.v(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
            if(field.Energy() == "on"){
                if (currentCell->type() == cell_type::ADIABATIC_WALL){
                    // T(i,j) = (T(i,j-1)+T(i-1,j))/2
                    field.setT(currentCell->i(),currentCell->j(),0.5*(field.T(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()) + field.T(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j())));            
                }
                if (currentCell->type() == cell_type::HOT_WALL || currentCell->type() == cell_type::COLD_WALL){
                    // T(i,j) = 2*T_wall - (T(i,j-1)+T(i-1,j))/2
                    field.setT(currentCell->i(),currentCell->j(),2*_wall_temperature - 0.5* (field.T(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()) + field.T(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j())));
                } 
            } 
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
    for(const auto& currentCell: _cells){
        if(currentCell->is_border(border_position::BOTTOM)){   //Top cells
            field.setv(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),0.0);
            // u doesn't lie on boundary, setting the avg value with neighbour (bottom) cell as the specified wall velocity
            field.setu(currentCell->i(),currentCell->j(),2*_wall_velocity.at(currentCell->wall_id())-(field.u(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j())));
            field.setp(currentCell->i(),currentCell->j(),field.p(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
            field.setg(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),field.v(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
        }
    }
}

InflowBoundary::InflowBoundary(std::vector<Cell *> cells, double inlet_velocity_x, double inlet_velocity_y) 
    : _cells(cells), _inlet_velocity_x(inlet_velocity_x), _inlet_velocity_y(inlet_velocity_y) {};

void InflowBoundary::apply(Fields &field) {
    for(const auto& currentCell: _cells){
        if(currentCell->is_border(border_position::RIGHT)){   //Fluid on right cells
            // U(i,j) = U_in
            field.setu(currentCell->i(),currentCell->j(), _inlet_velocity_x);
            // V(i,j) = 2*V_in - V(i+1,j)
            field.setv(currentCell->i(),currentCell->j(), 2*_inlet_velocity_y -(field.v(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j())));
            // P(i,j) = P(i+1,j)
            field.setp(currentCell->i(),currentCell->j(),field.p(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j()));
            // F(i,j) = U(i,j)
            field.setf(currentCell->i(),currentCell->j(),field.u(currentCell->i(),currentCell->j()));
        }
        if(currentCell->is_border(border_position::LEFT)){   //Fluid on left cells
            // U(i-1,j) = U_in
            field.setu(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),_inlet_velocity_x);
            // V(i,j) = 2*V_in - V(i-1,j)
            field.setv(currentCell->i(),currentCell->j(),2*_inlet_velocity_y -(field.v(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j())));
            // P(i,j) = P(i-1,j) 
            field.setp(currentCell->i(),currentCell->j(),field.p(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));
            // F(i-1,j) = U(i-1,j)
            field.setf(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),field.u(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));
        }
        if(currentCell->is_border(border_position::TOP)){   //Fluid on top cells
            // U(i,j) = 2*U_in - U(i,j+1)
            field.setu(currentCell->i(),currentCell->j(),2*_inlet_velocity_x - field.u(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()));
            // V(i,j) = V_in
            field.setv(currentCell->i(),currentCell->j(),_inlet_velocity_y);
            // P(i,j) = P(i,j+1)
            field.setp(currentCell->i(),currentCell->j(),field.p(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()));
            // G(i,j) = V(i,j)
            field.setg(currentCell->i(),currentCell->j(),field.v(currentCell->i(),currentCell->j()));
        }
        if(currentCell->is_border(border_position::BOTTOM)){   //Fluid on bottom cells
            // U(i,j) = 2*U_in - U(i,j-1)
            field.setu(currentCell->i(),currentCell->j(),2*_inlet_velocity_x -(field.u(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j())));
             // V(i,j-1) = V_in
            field.setv(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),_inlet_velocity_y);
            // P(i,j) = P(i,j-1)
            field.setp(currentCell->i(),currentCell->j(),field.p(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
            // G(i,j-1) = V(i,j-1)
            field.setg(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),field.v(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
        }
    }
}

OutflowBoundary::OutflowBoundary(std::vector<Cell *> cells, double outflow_pressure) 
    : _cells(cells), _outflow_pressure(outflow_pressure) {
    };

void OutflowBoundary::apply(Fields &field) {
    for(const auto& currentCell: _cells){
        if(currentCell->is_border(border_position::LEFT)){   //Fluid on left cells
            // U(i-1,j) = U(i-2,j) 
            field.setu(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),field.u(currentCell->neighbour(border_position::LEFT)->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->neighbour(border_position::LEFT)->j()));
            // V(i,j) = V(i-1,j)
            field.setv(currentCell->i(),currentCell->j(),field.v(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));
            // P(i,j) = 2*P_out - P(i-1,j) 
            field.setp(currentCell->i(),currentCell->j(),2 * _outflow_pressure - field.p(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));
            // F(i-1,j) = U(i-1,j)
            field.setf(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j(),field.u(currentCell->neighbour(border_position::LEFT)->i(),currentCell->neighbour(border_position::LEFT)->j()));
        }
        if(currentCell->is_border(border_position::RIGHT)){   //Fluid on right cells
            // U(i,j) = U(i+1,j) 
            field.setu(currentCell->i(),currentCell->j(),field.u(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j()));
            // V(i,j) = V(i+1,j)
            field.setv(currentCell->i(),currentCell->j(),field.v(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j()));
            // P(i,j) = 2*P_out - P(i+1,j)
            field.setp(currentCell->i(),currentCell->j(),2 * _outflow_pressure - field.p(currentCell->neighbour(border_position::RIGHT)->i(),currentCell->neighbour(border_position::RIGHT)->j()));
            // F(i,j) = U(i,j)
            field.setf(currentCell->i(),currentCell->j(),field.u(currentCell->i(),currentCell->j()));
        }
        if(currentCell->is_border(border_position::BOTTOM)){   //Fluid on bottom cells
            // V(i,j-1) = V(i,j-2)
            field.setv(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),field.v(currentCell->neighbour(border_position::BOTTOM)->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->neighbour(border_position::BOTTOM)->j()));
            // U(i,j) = U(i,j-1)
            field.setu(currentCell->i(),currentCell->j(),field.u(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
            // P(i,j) = 2*P_out - P(i,j-1)
            field.setp(currentCell->i(),currentCell->j(),2 * _outflow_pressure - field.p(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
            // G(i,j-1) = V(i,j-1)
            field.setg(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j(),field.v(currentCell->neighbour(border_position::BOTTOM)->i(),currentCell->neighbour(border_position::BOTTOM)->j()));
        }
        if(currentCell->is_border(border_position::TOP)){   //Fluid on top cells
            // V(i,j) = V(i,j+1)
            field.setv(currentCell->i(),currentCell->j(),field.v(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()));
            // U(i,j) = U(i,j+1)
            field.setu(currentCell->i(),currentCell->j(),field.u(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()));
            // P(i,j) = 2*P_out - P(i,j+1)
            field.setp(currentCell->i(),currentCell->j(),2 * _outflow_pressure - field.p(currentCell->neighbour(border_position::TOP)->i(),currentCell->neighbour(border_position::TOP)->j()));
            // G(i,j) = V(i,j)
            field.setg(currentCell->i(),currentCell->j(),field.v(currentCell->i(),currentCell->j()));
        }
    }
}
