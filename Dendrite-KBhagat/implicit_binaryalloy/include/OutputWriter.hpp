#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
class OutputWriter{
    public:
        explicit OutputWriter(const std::string& output_dir);

        //field_names: one name per component of the FESystem, e.g. {"u", "phi"}
        void write_vtu_and_pvd(const DoFHandler<Nsd>& dof_handler, const Vector<double>& solution,
                               const std::vector<std::string>& field_names,
                               const unsigned int timestep, std::ofstream &pvd_output);

        void initiate_pvd(std::ofstream &pvd_output);
        void finish_pvd(std::ofstream &pvd_output);

        const std::string& get_output_dir() const { return output_dir_; }//return output directory location

    private:
        std::string output_dir_;
};

#include "OutputWriter.tpp"
