#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
OutputWriter<Nsd,BfOrder>::OutputWriter(const std::string& output_dir)
: output_dir_(output_dir)
{
    std::filesystem::create_directories(output_dir_);
}

template <unsigned int Nsd, unsigned int BfOrder>
void OutputWriter<Nsd,BfOrder>::write_vtu_and_pvd(const DoFHandler<Nsd>& dof_handler, const Vector<double>& solution,
                                                  const std::vector<std::string>& field_names,
                                                  const unsigned int timestep, std::ofstream &pvd_output){
    DataOut<Nsd> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, field_names);

    data_out.build_patches();

    std::string filename = "solution-" + dealii::Utilities::int_to_string(timestep, 7) + ".vtu"; //7 is for padding integer names in the file names
    std::ofstream vtu_output(output_dir_ + "/" + filename);
    data_out.write_vtu(vtu_output);
    pvd_output << "    <DataSet timestep=\"" << timestep << "\" group=\"\" part=\"0\" file=\"" << filename << "\"/>\n"; //append to pvd file (same folder)
}

template <unsigned int Nsd, unsigned int BfOrder>
void OutputWriter<Nsd,BfOrder>::initiate_pvd(std::ofstream &pvd_output){
    pvd_output << "<?xml version=\"1.0\"?>\n"
               << "<VTKFile type=\"Collection\" version=\"0.1\" "
               << "byte_order=\"LittleEndian\" "
               << "compressor=\"vtkZLibDataCompressor\">\n"
               << "  <Collection>\n";
}

template <unsigned int Nsd, unsigned int BfOrder>
void OutputWriter<Nsd,BfOrder>::finish_pvd(std::ofstream &pvd_output){
    pvd_output << "  </Collection>\n" << "</VTKFile>\n";
    pvd_output.close();
}
