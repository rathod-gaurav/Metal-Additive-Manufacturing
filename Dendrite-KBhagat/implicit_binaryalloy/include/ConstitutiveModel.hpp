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

    private:
        double D_l_; double D_s_; double k_;
        double eps4_; double Gamma_g_; double m_l_; double c0_;
        double G_; double Vp_;
        double lambda0_; double a1_; double a2_;
};
