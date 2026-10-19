#pragma once

template<unsigned int Nsd, unsigned int Nne>
class PureMelt{
    public:
        PureMelt(
            const double x_ll, const double x_ul, const unsigned int n_elements
        );
        void run();
    
    private:
        const double x_ll_, x_ul_;
        const unsigned int n_elements_;

        void make_grid();

        Triangulation<Nsd> triangulation;
};

#include "PureMelt.tpp"
#include "MakeGrid.tpp"
