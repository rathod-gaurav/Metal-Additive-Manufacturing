#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
void PureMelt<Nsd,BfOrder>::make_grid(){
    GridGenerator::subdivided_hyper_cube(triangulation, n_elements_, x_ll_, x_ul_);
    std::cout << "Number of active cells: " << triangulation.n_active_cells() << std::endl;
}