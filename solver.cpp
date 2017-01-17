#include "solver.hpp"

#include <cassert>
#include <fstream>
#include <sstream>

Field::Field(size_t x_cells, size_t y_cells):
    base(x_cells * y_cells, 0),
    x_cells(x_cells),
    y_cells(y_cells) {}

void Field::resize(size_t x_cells, size_t y_cells) {
    
    this->x_cells = x_cells;
    this->y_cells = y_cells;
    
    base.resize(x_cells * y_cells);
}

double& Field::operator()(size_t x, size_t y) {
    
    assert((x < x_cells) && (y < y_cells));
    return base[x * y_cells + y];
}

double Field::operator()(size_t x, size_t y) const {
    
    assert((x < x_cells) && (y < y_cells));
    return base[x * y_cells + y];
}

Field Field::operator+(Field const &field) const {

    Field sum(*this);
    sum += field;

    return sum;
}

Field Field::operator-(Field const &field) const {

    Field diff(*this);
    diff -= field;

    return diff;
}

Field& Field::operator+=(Field const &field) {
    assert(x_cells == field.x_cells && y_cells == field.y_cells);

    size_t size = field.base.size();
    for (size_t i = 0; i < size; ++i) {
        base[i] += field.base[i];
    }

    return *this;
}

Field &Field::operator-=(Field const &field) {
    assert(x_cells == field.x_cells && y_cells == field.y_cells);

    size_t size = field.base.size();
    for (size_t i = 0; i < size; ++i) {
        base[i] -= field.base[i];
    }

    return *this;
}

Field Field::operator*(double multiplier) const {

    Field scaled(*this);
    for (double &element: scaled.base) {
        element *= multiplier;
    }

    return scaled;
}

const std::vector<double> &Field::getBase() const {
    return base;
}


ShallowWaterSolver::ShallowWaterSolver(
    size_t x_cells,
    size_t y_cells,
    double time_step,
    double dx,
    double gravity):

    dt(time_step),
    dx(dx),
    time_elapsed(0.0),
    x_cells(x_cells),
    y_cells(y_cells),
    water_height(x_cells, y_cells),
    vx(x_cells + 1, y_cells),
    vy(x_cells, y_cells + 1),
    surface_level(x_cells, y_cells),
    gravity(gravity) {}
    
void ShallowWaterSolver::run(size_t iterations) {
    
    for (size_t iteration = 0; iteration < iterations; ++iteration) {
        
//        rungeKutta4();
        euler();
        apply_reflecting_boundary_conditions();
//        output(iteration);
    }
}

Field ShallowWaterSolver::calculate_rhs_height() {
    
    Field new_height(x_cells, y_cells);

    for (size_t i = 1; i < x_cells - 1; ++i) {
        for (size_t j = 1; j < y_cells - 1; ++j) {
            double etha_loc = water_height(i, j);
            double vx_loc = (vx(i, j) + vx(i+1, j)) / 2;
            double vy_loc = (vy(i, j) + vy(i, j+1)) / 2;

            double dvxdx = (vx(i+1, j) - vx(i, j)) / dx;
            double dvydy = (vy(i, j+1) - vy(i, j)) / dx;

            double dethadx = (water_height(i+1, j) - water_height(i-1, j)) / (2 * dx);
            double dethady = (water_height(i, j+1) - water_height(i, j-1)) / (2 * dx);

            new_height(i, j) = -(dethadx * vx_loc + dethady * vy_loc + etha_loc * (dvxdx + dvydy));
        }
    }

    return new_height;
}

Field ShallowWaterSolver::calculate_rhs_vx() {
    
    // inner cells only
    Field new_vx(x_cells + 1, y_cells);

    for (size_t i = 1; i < x_cells; ++i) {
        for (size_t j = 1; j < y_cells - 1; ++j) {
            double vx_loc = vx(i, j);
            double vy_loc = (vy(i-1, j) + vy(i-1, j+1) + vy(i, j) + vy(i, j+1)) / 4;

            double dvxdx = (vx(i+1, j) - vx(i-1, j)) / (2 * dx);
            double dvxdy = (vx(i, j+1) - vx(i, j-1)) / (2 * dx);

            double dhdx = (surface_level(i, j) + water_height(i, j) - surface_level(i-1, j) - water_height(i-1, j)) / dx;

            new_vx(i, j) = -(vx_loc * dvxdx + vy_loc * dvxdy + gravity * dhdx);
        }
    }

    return new_vx;
}
            
Field ShallowWaterSolver::calculate_rhs_vy() {

    // inner cells only
    Field new_vy(x_cells, y_cells + 1);

    for (size_t i = 1; i < x_cells - 1; ++i) {
        for (size_t j = 1; j < y_cells; ++j) {
            double vx_loc = (vx(i, j) + vx(i+1, j) + vx(i, j-1) + vx(i+1, j-1))/ 4;
            double vy_loc = vy(i, j);

            double dvydx = (vy(i+1, j) - vy(i-1, j)) / (2 * dx);
            double dvydy = (vy(i, j+1) - vy(i, j-1)) / (2 * dx);

            double dhdy = (surface_level(i, j) + water_height(i, j) - surface_level(i, j-1) - water_height(i, j-1)) / dx;

            new_vy(i, j) = -(vx_loc * dvydx + vy_loc * dvydy + gravity * dhdy);
        }
    }

    return new_vy;
}

Field ShallowWaterSolver::advect_height(double time_step) const{

    Field advected_height(water_height);
//    #pragma omp parallel for shared(x_cells, y_cells)
    for (size_t i = 1; i < x_cells - 1; ++i) {
        for (size_t j = 1; j < y_cells - 1; ++j) {

            double x = (i + 0.5) * dx;
            double y = (j + 0.5) * dx;

            double v_x = (vx(i, j) + vx(i+1, j)) / 2;
            double v_y = (vy(i, j) + vy(i, j+1)) / 2;

            x -= v_x * time_step;
            y -= v_y * time_step;

            //  reflect from the boundaries
            double x_boundary = x_cells * dx;
            double y_boundary = y_cells * dx;

            if (x < 0) {
                x = -x;
            }
            if (x > x_boundary) {
                x = -x + 2 * x_boundary;
            }
            if (y < 0) {
                y = -y;
            }
            if (y > y_boundary) {
                y = -y + 2 * y_boundary;
            }

            //  interpolate the advected value
            //  At first, compute the advected indices on the mesh
            double advected_i = x / dx - 0.5;
            double advected_j = y / dx - 0.5;

            size_t advected_i_floor = (size_t) advected_i;
            size_t advected_j_floor = (size_t) advected_j;

            double i_offset = advected_i - advected_i_floor;
            double j_offset = advected_j - advected_j_floor;

            assert((i_offset >= 0.0) && (j_offset >= 0.0) && (i_offset <= 1.0) && (j_offset <= 1.0));

            advected_height(i, j) = water_height(advected_i_floor, advected_j_floor) * (1.0 - i_offset) * (1.0 - j_offset)
                         + water_height(advected_i_floor + 1, advected_j_floor) * i_offset * (1.0 - j_offset)
                         + water_height(advected_i_floor, advected_j_floor + 1) * (1.0 - i_offset) * j_offset
                         + water_height(advected_i_floor + 1, advected_j_floor + 1) * i_offset * j_offset;
        }
    }

    return advected_height;
}

Field ShallowWaterSolver::advect_vx(double time_step) const {

    Field advected_vx(vx);
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t i = 1; i < x_cells; ++i) {
        for (size_t j = 1; j < y_cells - 1; ++j) {

            double x = i * dx;
            double y = (j + 0.5) * dx;

            double v_x = vx(i, j);
            double v_y = (vy(i, j) + vy(i, j+1) + vy(i-1, j) + vy(i-1, j+1)) / 4;

            x -= v_x * time_step;
            y -= v_y * time_step;

            //  reflect from the boundaries
            double x_boundary = x_cells * dx;
            double y_boundary = y_cells * dx;

            if (x < 0) {
                x = -x;
            }
            if (x > x_boundary) {
                x = -x + 2 * x_boundary;
            }
            if (y < 0) {
                y = -y;
            }
            if (y > y_boundary) {
                y = -y + 2 * y_boundary;
            }

            //  interpolate the advected value
            //  At first, compute the advected indices on the mesh
            double advected_i = x / dx;
            double advected_j = y / dx - 0.5;

            size_t advected_i_floor = (size_t) advected_i;
            size_t advected_j_floor = (size_t) advected_j;

            double i_offset = advected_i - advected_i_floor;
            double j_offset = advected_j - advected_j_floor;

            assert((i_offset >= 0.0) && (j_offset >= 0.0) && (i_offset <= 1.0) && (j_offset <= 1.0));

            advected_vx(i, j) = vx(advected_i_floor, advected_j_floor) * (1.0 - i_offset) * (1.0 - j_offset)
                         + vx(advected_i_floor + 1, advected_j_floor) * i_offset * (1.0 - j_offset)
                         + vx(advected_i_floor, advected_j_floor + 1) * (1.0 - i_offset) * j_offset
                         + vx(advected_i_floor + 1, advected_j_floor + 1) * i_offset * j_offset;
        }
    }

    return advected_vx;
}

Field ShallowWaterSolver::advect_vy(double time_step) const {

    Field advected_vy(vy);
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t i = 1; i < x_cells - 1; ++i) {
        for (size_t j = 1; j < y_cells; ++j) {

            double x = (i + 0.5) * dx;
            double y = j * dx;

            double v_x = (vx(i, j) + vx(i+1, j) + vx(i, j-1) + vx(i+1, j-1)) / 4;
            double v_y = vy(i, j);

            x -= v_x * time_step;
            y -= v_y * time_step;

            //  reflect from the boundaries
            double x_boundary = x_cells * dx;
            double y_boundary = y_cells * dx;

            if (x < 0) {
                x = -x;
            }
            if (x > x_boundary) {
                x = -x + 2 * x_boundary;
            }
            if (y < 0) {
                y = -y;
            }
            if (y > y_boundary) {
                y = -y + 2 * y_boundary;
            }

            //  interpolate the advected value
            //  At first, compute the advected indices on the mesh
            double advected_i = x / dx - 0.5;
            double advected_j = y / dx;

            size_t advected_i_floor = (size_t) advected_i;
            size_t advected_j_floor = (size_t) advected_j;

            double i_offset = advected_i - advected_i_floor;
            double j_offset = advected_j - advected_j_floor;

            assert((i_offset >= 0.0) && (j_offset >= 0.0) && (i_offset <= 1.0) && (j_offset <= 1.0));

            advected_vy(i, j) = vy(advected_i_floor, advected_j_floor) * (1.0 - i_offset) * (1.0 - j_offset)
                         + vy(advected_i_floor + 1, advected_j_floor) * i_offset * (1.0 - j_offset)
                         + vy(advected_i_floor, advected_j_floor + 1) * (1.0 - i_offset) * j_offset
                         + vy(advected_i_floor + 1, advected_j_floor + 1) * i_offset * j_offset;
        }
    }

    return advected_vy;
}

void ShallowWaterSolver::update_height(double time_step) {

//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t i = 1; i < x_cells - 1; ++i) {
        for (size_t j = 1; j < y_cells - 1; ++j) {
            
            // divergence from the staggered grid
            double divergence = (vx(i+1,j) - vx(i,j)
                                 + vy(i,j+1) - vy(i,j)) / dx;
            water_height(i,j) -= water_height(i,j) * divergence * time_step;
        }
    }
}

void ShallowWaterSolver::update_vx(double time_step) {

//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t j = 1; j < y_cells - 1; ++j) {
        for (size_t i = 2; i < x_cells - 1; ++i) {
            
            double height_above_zero_left = surface_level(i - 1, j) + water_height(i - 1, j);
            double height_above_zero_right = surface_level(i, j) + water_height(i, j);
            
            vx(i,j) += gravity * (height_above_zero_left - height_above_zero_right) / dx * time_step;
        }
    }
}

void ShallowWaterSolver::update_vy(double time_step) {
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t j = 2; j < y_cells - 1; ++j) {
        for (size_t i = 1; i < x_cells - 1; ++i) {

            double height_above_zero_down = surface_level(i, j - 1) + water_height(i, j - 1);
            double height_above_zero_up = surface_level(i, j) + water_height(i, j);

            vy(i,j) += gravity * (height_above_zero_down - height_above_zero_up) / dx * time_step;
        }
    }
}

void ShallowWaterSolver::apply_reflecting_boundary_conditions() {
    
    //  left boundary
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t j = 0; j < y_cells; ++j) {
        water_height(0, j) = water_height(1, j) + surface_level(1, j) - surface_level(0, j);
        
        vx(1, j) = 0.0;
        vy(0, j) = 0.0;
//        vx(0, j) = 0.0; // never accessed
    }
    
    //  down boundary
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t i = 0; i < x_cells; ++i) {
        water_height(i, 0) = water_height(i, 1);
        
        vy(i, 1) = 0.0;
        vx(i, 0) = 0.0;
//        vy(i, 0) = 0.0;
    }
    
    //  right boundary
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t j = 0; j < y_cells; ++j) {
        water_height(x_cells - 1, j) = water_height(x_cells - 2, j) + surface_level(x_cells - 2, j) - surface_level(x_cells - 1, j);
        
        vx(x_cells - 1, j) = 0.0;
//        vx(x_cells, j) = 0.0;
        vy(x_cells - 1, j) = 0.0;
    }
    
    //  upper boundary
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t i = 0; i < x_cells; ++i) {
        water_height(i, y_cells - 1) = water_height(i, y_cells - 2) + surface_level(i, y_cells - 2) - surface_level(i, y_cells - 1);
        
        vy(i, y_cells - 1) = 0.0;
        vx(i, y_cells - 1) = 0.0;
//        vy(i, y_cells) = 0.0;
    }
}

double ShallowWaterSolver::total_mass() {

    double mass = 0;
    for (size_t i = 0; i < x_cells; ++i) {
        for (size_t j = 0; j < y_cells; ++j) {
            mass += water_height(i, j);
        }
    }
    return mass;
}

double ShallowWaterSolver::total_energy() {

    double kinetic_energy = 0.0;
    double potential_energy = 0.0;

    for (size_t i = 0; i < x_cells; ++i) {
        for (size_t j = 0; j < y_cells; ++j) {
            double vx_loc = (vx(i,j) + vx(i+1,j)) / 2;
            double vy_loc = (vy(i,j) + vy(i,j+1)) / 2;
            double h = water_height(i,j);

            kinetic_energy += h * (vx_loc * vx_loc + vy_loc * vy_loc);
            potential_energy += h * h / 2;
        }
    }

    return kinetic_energy + potential_energy * gravity;
}


void ShallowWaterSolver::initialize_water_height(const std::vector<double> &input) {
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t i = 0; i < x_cells; ++i) {
        for (size_t j = 0; j < y_cells; ++j) {
            water_height(i, j) = input[i * x_cells + j];
        }
    }
}

void ShallowWaterSolver::initialize_vx(const std::vector<double> &input) {
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t i = 0; i <= x_cells; ++i) {
        for (size_t j = 0; j < y_cells; ++j) {
            vx(i, j) = input[i * x_cells + j];
        }
    }
}

void ShallowWaterSolver::initialize_vy(const std::vector<double> &input) {
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t i = 0; i < x_cells; ++i) {
        for (size_t j = 0; j <= y_cells; ++j) {
            vy(i, j) = input[i * x_cells + j];
        }
    }
}

void ShallowWaterSolver::initialize_surface_level(const std::vector<double> &input) {
//#pragma omp parallel for shared(x_cells, y_cells)
    for (size_t i = 0; i < x_cells; ++i) {
        for (size_t j = 0; j < y_cells; ++j) {
            surface_level(i, j) = input[i * x_cells + j];
        }
    }
}

void ShallowWaterSolver::rungeKutta4() {

    Field original_height(water_height);
    Field original_vx(vx);
    Field original_vy(vy);

    Field k1_height(advect_height(dt / 2.0));
    Field k1_vx(advect_vx(dt / 2.0));
    Field k1_vy(advect_vy(dt / 2.0));

    water_height = k1_height;
    vx = k1_vx;
    vy = k1_vy;

    update_height(dt / 2.0);
    update_vx(dt / 2.0);
    update_vy(dt / 2.0);

    k1_height = water_height - original_height;
    k1_vx = vx - original_vx;
    k1_vy = vy - original_vy;
    // k1 = h/2 * k1

    Field k2_height(advect_height(dt / 2.0));
    Field k2_vx(advect_vx(dt / 2.0));
    Field k2_vy(advect_vy(dt / 2.0));

    water_height = k2_height;
    vx = k2_vx;
    vy = k2_vy;

    update_height(dt / 2.0);
    update_vx(dt / 2.0);
    update_vy(dt / 2.0);

    // Leave only y = y_original + h/2 * k2
    water_height -= k1_height;
    vx -= k1_vx;
    vy -= k1_vy;

    k2_height = water_height - original_height;
    k2_vx = vx - original_vx;
    k2_vy = vy - original_vy;
    // k2 = h/2 * k2

    Field k3_height(advect_height(dt));
    Field k3_vx(advect_vx(dt));
    Field k3_vy(advect_vy(dt));

    water_height = k3_height;
    vx = k3_vx;
    vy = k3_vy;

    update_height(dt);
    update_vx(dt);
    update_vy(dt);
    // Now they are y = y_original + h/2 * k2 + h * k3

    water_height -= k2_height;
    vx -= k2_vx;
    vy -= k2_vy;

    k3_height = water_height - original_height;
    k3_vx = vx - original_vx;
    k3_vy = vy - original_vy;
    // k3 = h * k3

    Field k4_height(advect_height(dt / 2.0));
    Field k4_vx(advect_vx(dt / 2.0));
    Field k4_vy(advect_vy(dt / 2.0));

    water_height = k4_height;
    vx = k4_vx;
    vy = k4_vy;

    update_height(dt / 2.0);
    update_vx(dt / 2.0);
    update_vy(dt / 2.0);
    // Now they are y = y_original + h * k3 + h/2 * k4

    water_height -= k3_height;
    vx -= k3_vx;
    vy -= k3_vy;

    k4_height = water_height - original_height;
    k4_vx = vx - original_vx;
    k4_vy = vy - original_vy;
    // k4 = h/2 * k4

    water_height -= k4_height;
    vx -= k4_vx;
    vy -= k4_vy;

//    water_height += (k1_height*(1.0/3.0)  + k2_height*(2.0/3.0) + k3_height *(1.0/3.0) + k4_height*(2.0/6.0));
    water_height += ((k1_height + k3_height + k4_height)*(1.0/3.0)  + k2_height*(2.0/3.0));
    vx = ((k1_vx + k3_vx + k4_vx)*(1.0/3.0)  + k2_vx*(2.0/3.0) );
    vy = ((k1_vy + k3_vy + k4_vy)*(1.0/3.0)  + k2_vy*(2.0/3.0));

    time_elapsed += dt;
}

double ShallowWaterSolver::getTime_elapsed() const {
    return time_elapsed;
}

const Field &ShallowWaterSolver::getWater_height() const {
    return water_height;
}

const Field &ShallowWaterSolver::getSurface_level() const {
    return surface_level;
}

void ShallowWaterSolver::output(size_t iteration) {
    std::stringstream ss;
    ss << "output/height" << iteration << ".txt";
    std::ofstream fout(ss.str());
    for (size_t x = 0; x < x_cells; ++x) {
        for (size_t y = 0; y < y_cells; ++y) {
            fout << water_height(x, y) << ' ';
        }
        fout << '\n';
    }
}

void ShallowWaterSolver::euler() {
    Field new_height(advect_height(dt));
    Field new_vx(advect_vx(dt));
    Field new_vy(advect_vy(dt));

    water_height = new_height; //advect_height(dt);
    vx = new_vx;//advect_vx(dt);
    vy = new_vy;//advect_vy(dt);//

    update_height(dt);
    update_vx(dt);
    update_vy(dt);

    time_elapsed += dt;
}
