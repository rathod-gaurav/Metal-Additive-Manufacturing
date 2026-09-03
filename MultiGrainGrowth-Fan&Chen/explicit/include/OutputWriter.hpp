#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
class OutputWriter{
    public:
        explicit OutputWriter(const std::string& output_dir);

        void write_vtu(DoFHandler<Nsd>& dof_handler, Vector<double>& solution1, const unsigned int timestep);
        void write_pvd();

        const std::string& get_output_dir() const { return output_dir_; }//return output directory location

    private:
        std::string output_dir_;
};

#include "OutputWriter.tpp"