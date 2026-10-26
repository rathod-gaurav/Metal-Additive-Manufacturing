#pragma once

// QPFields is a container. for one element, it holds every field’s value and gradient at every quadrature point, at the new time level and at the old one.
// evaluate_fields() fills it, and element_residual() reads it.

template<unsigned int Nsd, typename Number>
struct QPFields{//everything the residual needs at the quadrature point of each element
    QPFields(const unsigned int n_fields, const unsigned int n_q)
    :   
        val(n_fields, std::vector<Number>(n_q)),
        grad(n_fields, std::vector<Tensor<1,Nsd,Number>>(n_q)), //Create a dynamic list of rank-1 tensors (each sized for Nsd dimensions using double precision), and size the list to hold n_q items  
        val_old(n_fields, std::vector<double>(n_q)),
        grad_old(n_fields, std::vector<Tensor<1,Nsd,double>>(n_q))
    {}

    std::vector<std::vector<Number>>                val;
    std::vector<std::vector<Tensor<1,Nsd,Number>>>  grad;
    std::vector<std::vector<double>>                val_old;
    std::vector<std::vector<Tensor<1,Nsd,double>>>  grad_old;
};