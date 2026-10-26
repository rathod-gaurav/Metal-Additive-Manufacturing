#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
template<typename Number>
void PureMelt<Nsd,BfOrder>::element_residual(
    const FEValues<Nsd>& fe_values, const QPFields<Nsd,Number>& fields, std::vector<Number>& Rlocal
) const {
    
}