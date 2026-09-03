#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
OutputWriter<Nsd,BfOrder>::OutputWriter(const std::string& output_dir)
: output_dir_(output_dir)
{
    std::filesystem::create_directories(output_dir_);
}

template <unsigned int Nsd, unsigned int BfOrder>
void OutputWriter<Nsd,BfOrder>::write_vtu(DoFHandler<Nsd>& dof_handler, Vector<double>& solution1, const unsigned int timestep){
    DataOut<Nsd> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution1, "phi");

    data_out.build_patches();

    std::string filename = output_dir_ + "/solution-" + dealii::Utilities::int_to_string(timestep, 7) + ".vtu"; //6 is for padding integer names in the file names
    std::ofstream vtu_output(filename);
    data_out.write_vtu(vtu_output);
}

template <unsigned int Nsd, unsigned int BfOrder>
void OutputWriter<Nsd,BfOrder>::write_pvd(){
    std::string pvd_filename = "final_output.pvd";
    std::ofstream pvd_output("final_solution.pvd");

    pvd_output << "<?xml version=\"1.0\"?>\n"
               << "<VTKFile type=\"Collection\" version=\"0.1\" "
               << "byte_order=\"LittleEndian\" "
               << "compressor=\"vtkZLibDataCompressor\">\n"
               << "  <Collection>\n";

    // Collect and sort all .vtu files in output_dir_ to ensure chronological order
    std::vector<fs::path> vtu_files;
    if (fs::exists(output_dir_) && fs::is_directory(output_dir_)) {
        for (const auto& entry : fs::directory_iterator(output_dir_)) {
            if (entry.path().extension() == ".vtu") {
                vtu_files.push_back(entry.path());
            }
        }
    }
    // std::sort(vtu_files.begin(), vtu_files.end());

    // Loop through each discovered file and append it to the PVD collection
    unsigned int timestep = 0;
    for (size_t i = 0; i < vtu_files.size(); ++i) {
        std::string filename = vtu_files[i].filename().string();
        pvd_output << "    <DataSet timestep=\"" << timestep << "\" group=\"\" part=\"0\" file=\"" << filename << "\"/>\n";
        timestep += 1; 
    }

    pvd_output << "  </Collection>\n" << "</VTKFile>\n";
    pvd_output.close();

    std::cout << "Output file written to " << pvd_filename << "." << std::endl;
}

