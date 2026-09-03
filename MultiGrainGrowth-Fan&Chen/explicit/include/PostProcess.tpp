#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::post_process(){
    //Assemble function phi as explained in the paper
    for (unsigned int k = 0; k < p_; ++k){
        const double *const row = &eta_np1[k][0];
        for (unsigned int j = 0; j < Nt; ++j){
            phi(j) += row[j] * row[j];
        }
    }
}