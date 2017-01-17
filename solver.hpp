#pragma once

// Finite difference approach
// 

//  Velocity and pressure location on a staggered grid:
//  +-- vy-- + -- vy-- +
//  |        |         |
//  |   h    vx   h    vx
//  |        |         |
//  +------- + ------- +
//
//  This enables us to prevent instabilities while computing symmetric spatial gradients.
//  Note that the gradients of velocities affect only pressure and vise versa.

#include <vector>
#include <stddef.h>
//#include <fstream>
//#include <sstream>
#include <omp.h>

    
class Field {
public:
    
    Field(size_t x_cells = 0, size_t y_cells = 0);
    Field(Field const & copy_from) = default;
    Field(Field&& move_from) = default;

    double &operator()(size_t x, size_t y);
    double  operator()(size_t x, size_t y) const;
    
    void resize(size_t x_cells, size_t y_cells);

    Field& operator=(Field const &copy_from) = default;

    Field operator+(Field const &field) const;
    Field operator-(Field const &field) const;
    Field& operator+=(Field const &field);
    Field& operator-=(Field const &field);
    Field operator*(double multiplier) const;

    const std::vector<double> &getBase() const;

private:
    
    std::vector<double> base;

    size_t x_cells;
    size_t y_cells;
};


//  Solution of the shallow water equation consists of 7 steps
//  1.  Advect height above zero level
//  2.  Advect v_x
//  3.  Advect v_y
//  4.  Update height (rhs)
//  5.  Compute h = eta' + g, height above surface level
//  6.  Update velocities
//  7.  Apply boundary conditions

class ShallowWaterSolver {
public:
    
    ShallowWaterSolver(size_t x_cells, size_t y_cells, double time_step, double dx, double gravity = 9.8);

    void initialize_water_height(const std::vector<double> &input);
    void initialize_vx(const std::vector<double> &input);
    void initialize_vy(const std::vector<double> &input);
    void initialize_surface_level(const std::vector<double> &input);
    
    void run(size_t iterations);

    const Field &getWater_height() const;
    const Field &getSurface_level() const;

    double getTime_elapsed() const;

    void output(size_t iteration);

    double total_mass();
    double total_energy();
private:

    void rungeKutta4();
    void euler();

//  Simple group of methods
    Field calculate_rhs_height();
    Field calculate_rhs_vx();
    Field calculate_rhs_vy();

//  Advanced group of methods
    Field advect_height(double time_step) const;
    Field advect_vx(double time_step) const;
    Field advect_vy(double time_step) const;
    void update_height(double time_step);
    void update_vx(double time_step);
    void update_vy(double time_step);

    void apply_reflecting_boundary_conditions();
    
    double dt;
    double dx;
    double time_elapsed;

    size_t x_cells;
    size_t y_cells;

    // Water cylinder height above the surface level
    Field water_height;
    Field vx;
    Field vy;

    Field surface_level;
    
    double const gravity;
};
