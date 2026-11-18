#pragma once

template<unsigned int Nsd>
struct NucleationSeed{
    Point<Nsd> center;
    double radius;
};

template <unsigned int Nsd>
std::vector<NucleationSeed<Nsd>> generate_random_seeds(const unsigned int Nseeds, const double x_ll, const double x_ul)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(x_ll, x_ul);
    std::vector<NucleationSeed<Nsd>> seeds(Nseeds);
    for (unsigned int i = 0; i < Nseeds; ++i)
    {
        if constexpr(Nsd==2){
            seeds[i].center = Point<Nsd>(dist(gen), dist(gen));
        }
        else if constexpr(Nsd==3){
            seeds[i].center = Point<Nsd>(dist(gen), dist(gen), dist(gen));
        }
        seeds[i].radius  = 5.0; //same radius for all seeds
    }
    return seeds;
}

template<unsigned int Nsd>
class InitialConditions : public Function<Nsd>{ //this class is the function in dealii's sense
    public: //everything public in function stays public in our class
        //class constructor
        InitialConditions(const double supersaturation,
                          const std::vector<NucleationSeed<Nsd>>& seeds,
                          const double lambda0 = 1.0)//we fix lambda0 specifically for the initial condition seed here
        : Function<Nsd>(2), //2 for two fields - temp u and phi
          supersaturation_(supersaturation),
          seeds_(seeds),
          lambda0_(lambda0)
        {}

        double value(const Point<Nsd>& p, const unsigned int component) const override{
            if(component == 0){ //component 0 is composition uc
                return supersaturation_; //initialize composition with supersaturation
            }
            // For component 1 (phase field phi), find the distance to the nearest seed
            double min_r = 1e20;
            double active_radius = 5.0;

            for (const auto& seed : seeds_) {
                double r = p.distance(seed.center);
                if (r < min_r) {
                    min_r = r;
                    active_radius = seed.radius;
                }
            }
            return -std::tanh((min_r - active_radius)/(std::sqrt(2.0)*lambda0_)); //initialize phase field with a tanh profile of radius r about seed_center_
        }

    private:
        const double supersaturation_; std::vector<NucleationSeed<Nsd>> seeds_; const double lambda0_;
};