#pragma once

class ConstitutiveModel{
    public:
        ConstitutiveModel(const double D, const double eps4, const double lambda0,
                          const double tau0, const double xi)
        : D_(D), eps4_(eps4), lambda0_(lambda0), tau0_(tau0), xi_(xi) {}

        double diffusivity() const { return D_; }

        //bulk driving force  df/dphi = g'(phi) + xi u p'(phi) = -phi + phi^3 + xi u (1 - phi^2)^2
        template <typename T>
        T bulk_driving_force(const T& phi, const T& u) const{
            const T one_m_phi2 = 1.0 - phi*phi;
            return -phi + phi*phi*phi + xi_*u*one_m_phi2*one_m_phi2;
        }

        //latent heat function h(phi) and its slope h'(phi) (U = u - h(phi)/2); here h = phi
        template <typename T>
        T latent_heat_slope(const T& /*phi*/) const{ return T(1.0); }

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

        //interface coefficients from the anisotropy:
        //lambda(n) = lambda0 a_s,  d lambda / d q,  tau(n) = tau0 a_s^2
        template <int dim, typename T>
        void interface_coefficients(const Tensor<1,dim,T>& q, T& lambda, Tensor<1,dim,T>& dlambda_dq, T& tau) const{
            T a;
            Tensor<1,dim,T> da_dq;
            anisotropy(q, a, da_dq);
            lambda     = lambda0_*a;
            dlambda_dq = lambda0_*da_dq;
            tau        = tau0_*a*a;
        }

    private:
        const double D_, eps4_, lambda0_, tau0_, xi_;
};
