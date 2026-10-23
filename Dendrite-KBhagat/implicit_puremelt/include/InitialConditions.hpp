#pragma once

template<unsigned int Nsd>
class InitialConditions : public Function<Nsd>{ //this class is the function in dealii's sense
    public: //everything public in function stays public in our class
        //class constructor
        InitialConditions(const double undercooling,
                          Point<Nsd>& seed_center,
                          const double seed_radius,
                          const double lambda0 = 1.0)//we fix lambda0 specifically for the initial condition seed here
        : Function<Nsd>(2), //2 for two fields - temp u and phi
          undercooling_(undercooling),
          seed_center_(seed_center),                    
          seed_radius_(seed_radius),
          lambda0_(lambda0)
        {}

        double value(const Point<Nsd>& p, const unsigned int component) const override{
            if(component == 0){ //component 0 is temperature u
                return undercooling_; 
            }
            const double r = p.distance(seed_center_);        
            return -std::tanh((r - seed_radius_)/(std::sqrt(2.0)*lambda0_)); //initialize phase field with a tanh profile of radius r about seed_center_
        }

    private:
        const double undercooling_; const Point<Nsd> seed_center_; const double seed_radius_, lambda0_;
};