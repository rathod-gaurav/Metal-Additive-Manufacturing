#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::make_grid(){
    GridGenerator::hyper_cube(triangulation, x_ll_, x_ul_, true); //true means colorize=true, It colors the boundary faces with boundary_id = 0,1,2,3
    triangulation.refine_global(n_refinements_);

    std::cout << "Number of active cells: " << triangulation.n_active_cells() << std::endl;
}