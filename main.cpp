#include "solver.hpp"
#include "render.hpp"
#include <fstream>
#include <cassert>
#include <vector>
#include <string>

void load_field_data(std::string const file_name, std::vector<double> &field, size_t &x_cells, size_t &y_cells) {
    /*
        Data file format follows:

        x_cells y_cells
        f_00 f_01 f_02 ... f_{0 y_cells-1}
        ...
        f_{x_cells-1 0} ... f_{x_cells-1 y_cells-1}
    */

    std::ifstream input("input/" + file_name);
    input.exceptions(std::istream::failbit | std::istream::badbit);

    input >> x_cells >> y_cells;

    field.resize(x_cells * y_cells);
    for (double &value: field) {
        input >> value;
    }
}

int main() {

    size_t x_cells;
    size_t y_cells;

    std::vector<double> initial_height;
    std::vector<double> initial_velocity_x;
    std::vector<double> initial_velocity_y;
    std::vector<double> initial_surface_level;

    load_field_data("velocity_x.txt", initial_velocity_x, x_cells, y_cells);
    load_field_data("velocity_y.txt", initial_velocity_y, x_cells, y_cells);
    load_field_data("height.txt", initial_height, x_cells, y_cells);
    load_field_data("surface.txt", initial_surface_level, x_cells, y_cells);
    assert(x_cells * y_cells == initial_height.size());
    // Note that now x_cells and y_cells represent the central field values -- exactly what is needed

    // Courant scheme stability requirement is time_step < dx (well...)
    double const time_step = 0.01;
    double const dx = 0.1;

    ShallowWaterSolver solver(x_cells, y_cells, time_step, dx);

    solver.initialize_water_height(initial_height);
    solver.initialize_vx(initial_velocity_x);
    solver.initialize_vy(initial_velocity_y);
    solver.initialize_surface_level(initial_surface_level);

    std::vector<GLfloat> water_height(initial_height.begin(), initial_height.end());
    std::vector<GLfloat> surface_level(initial_surface_level.begin(), initial_surface_level.end());
    for (size_t i = 0; i < water_height.size(); ++i) {
        water_height[i] += surface_level[i];
    }

    VisualEngine visual_engine(x_cells, y_cells, dx, 1200, 900);
    visual_engine.update_vertex_values(&water_height, &surface_level);

//    solver.run(1600);
    while(!visual_engine.should_stop()) {

//    for (size_t iteration = 0; iteration < 600; ++iteration) {
//#pragma omp parallel for
//        for (int threadId = 0; threadId < 2; ++threadId) {
//            if (threadId == 0) {
                visual_engine.render();
//            } else {
//                solver.run(1);
//            }
//        }
        if (visual_engine.start_simulation()) {
            std::cout << solver.total_energy() << '\n';
            solver.run(1);
            std::vector<double> const *wh = &solver.getWater_height().getBase();
            water_height = std::vector<GLfloat>(wh->begin(), wh->end());
            for (size_t i = 0; i < water_height.size(); ++i) {
                water_height[i] += surface_level[i];
            }
            visual_engine.update_vertex_values(&water_height);
        }
//        solver.copy_height_to(visual_engine.height);
//        visual_engine.render();
    }

    return 0;
}
