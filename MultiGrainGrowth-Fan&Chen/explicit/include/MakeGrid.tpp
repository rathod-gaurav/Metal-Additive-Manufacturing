#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::make_grid(){
    GridGenerator::hyper_cube(triangulation, x_ll_, x_ul_);
    triangulation.refine_global(n_refinements_);

    std::cout << "Number of active cells: " << triangulation.n_active_cells() << std::endl;
}