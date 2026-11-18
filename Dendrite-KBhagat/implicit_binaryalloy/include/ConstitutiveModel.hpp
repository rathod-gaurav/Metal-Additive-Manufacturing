#pragma once

class ConstitutiveModel{
    public:
        ConstitutiveModel(
            const double D_l, const double D_s, const double k,
            const double eps4, const double Gamma_g, const double m_l, const double c0,
            const double G, const double Vp,
            const double lambda0, const double a1, const double a2
        ) : D_l_(D_l), D_s_(D_s), k_(k), 
            eps4_(eps4), Gamma_g_(Gamma_g), m_l_(m_l), c0_(c0), 
            G_(G), Vp_(Vp), 
            lambda0_(lambda0), a1_(a1), a2_(a2)
        {}

        double xi() const {return a1_;}
        double tau0() const {return (a1_*a2_*xi()*lambda0_*lambda0_)/D_l_;}
        double diffusivity_non_dim() const {return (D_l_*tau0())/(lambda0_*lambda0_);}
        double partition_coefficient() const {return k_;}
        double lambda0() const {return lambda0_;}

        //four-fold anisotropy a_s(q) = (1-3eps4) + 4 eps4 (qx^4+qy^4)/|q|^4, q = grad phi
        //returns a_s and d a_s / d q
        template <int dim, typename T>
        void anisotropy(const Tensor<1,dim,T>& q, T& a, Tensor<1,dim,T>& da_dq) const{
            static_assert(dim == 2, "anisotropy() is written for 2D");
            const double delta = 1e-12; //regularisation where grad phi = 0 (bulk)
            const T qx = q[0], qy = q[1];
            const T q2 = qx*qx + qy*qy + delta;
            const T q4 = q2*q2;
            const T q6 = q4*q2;
            const T S  = (qx*qx*qx*qx + qy*qy*qy*qy)/q4;
            a = (1.0 - 3.0*eps4_) + 4.0*eps4_*S;
            da_dq[0] = 4.0*eps4_*(4.0*qx*qy*qy*(qx*qx - qy*qy)/q6);
            da_dq[1] = 4.0*eps4_*(4.0*qx*qx*qy*(qy*qy - qx*qx)/q6);
        }

    private:
        double D_l_; double D_s_; double k_;
        double eps4_; double Gamma_g_; double m_l_; double c0_;
        double G_; double Vp_;
        double lambda0_; double a1_; double a2_;
};
