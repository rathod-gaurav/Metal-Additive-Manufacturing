#pragma once
//function to apply periodic boundary conditions on the global system matrices and vectors

template <unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::apply_bc(){
    //Stage1: The triangulation
    std::vector<GridTools::PeriodicFacePair<typename Triangulation<Nsd>::cell_iterator>> matched_pairs; //An empty container to be filled. A PeriodicFacePair records one matched pair: which cell and face number on side A, which on side B, and the orientation relating them.
    GridTools::collect_periodic_faces(triangulation, 0, 1, 0, matched_pairs); //collect periodic faces between boundary_id 0 and 1 in direction 0 (x-direction)
    GridTools::collect_periodic_faces(triangulation, 2, 3, 1, matched_pairs); //collect periodic faces between boundary_id 2 and 3 in direction 1 (y-direction)
    // after both calls matched_pairs holds every match in both directions.

    triangulation.add_periodicity(matched_pairs); //add periodicity to the triangulation //Now the mesh itself knows cells on opposite edges are neighbours. This affects refinement (a cell refined on one edge forces its periodic partner to match) and, in parallel, ensures periodic partners land on communicating processes.

    //Stage2: The constraints
    // The identical calls, but on the DoFHandler and producing DoFHandler-typed iterators. The geometric matching is the same; what differs is that these iterators can be asked which DoFs live on a face. Stage 1's pairs can't do that, which is why you can't reuse them.
    std::vector<GridTools::PeriodicFacePair<typename DoFHandler<Nsd>::cell_iterator>> dof_pairs;
    GridTools::collect_periodic_faces(dof_handler, 0, 1, 0, dof_pairs);
    GridTools::collect_periodic_faces(dof_handler, 2, 3, 1, dof_pairs);
    constraints.clear();
    DoFTools::make_periodicity_constraints<int(Nsd),int(Nsd)>(dof_pairs, constraints); //The actual work. For each matched face pair, walk the DoFs on both faces and add a constraint saying "DoF on side B equals DoF on side A." After this, the B-side DoFs are constrained — they're no longer independent unknowns.
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); //Add hanging-node constraints into the same object (they must share one AffineConstraints, since a DoF can be both periodic and hanging)
    constraints.close(); //sorts and resolves chains — a DoF constrained to a DoF that is itself constrained gets flattened. Nothing works until close() is called.
    
}