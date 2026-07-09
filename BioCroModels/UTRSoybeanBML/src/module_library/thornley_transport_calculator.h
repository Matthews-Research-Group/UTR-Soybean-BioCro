#ifndef THORNLEY_TRANSPORT_CALCULATOR_H
#define THORNLEY_TRANSPORT_CALCULATOR_H

#include <vector>
#include <string>
#include "../framework/module.h"
#include "../framework/state_map.h"
#include "thornley_nutrient_dynamics.h"

/**
 * @class thornley_transport_calculator
 * 
 * @brief Part of the Thornley model for plant growth; although these calculations
 * could be calculated by the multi_organ_transport derivative module in principle, we
 * perform them here in a steady state model so they are available as outputs of a
 * BioCro simulation.
 * 
 * This model and its associated modules are primarily based on Thornley, J. H. M.
 * "A Model to Describe the Partitioning of Photosynthate during Vegetative Plant
 * Growth." Ann Bot 36, 419–430 (1972).
 * 
 * In this model, two types of carbon are distinguished (`substrate` carbon and
 * `structural` carbon) and each organ has its own separate pools of substrate and
 * structural carbon.
 * 
 * Substrate carbon is ultimately derived from photosynthesis. It is utilized for
 * growth and respiration, and it is transported between organs.
 * 
 * Structural carbon is lost during senescene, with some of the lost mass going to
 * the litter and the remainder being reclaimed as substrate carbon. Moreover, some
 * organs (called `storage organs`) can convert structural carbon to substrate carbon
 * in a process that is distinct from senescence.
 * 
 * A key assumption of this model is that substrate is transported pairwise between
 * organs with a rate determined by the difference in substrate concentration between
 * those two organs, i.e.,
 * 
 * transport_rate = conductance * (difference in substrate concentration).
 * 
 * As before, `mass_fraction` is used as a replacement for substrate concentration. Thornley
 * also points out that it may be necessary to rescale the conductance as the plant grows and
 * suggests using a scaling factor `beta` given by a power of the total amount of structural
 * carbon in the plant, i.e.,
 * 
 * beta = (total structural carbon)^m,
 * 
 * where `m` is the "beta exponent" and is typically either 1.0 (which results in "steady state
 * exponential growth") or 0.33 (corresponding to "hydraulic resistance"). Thus, the final equation
 * for the transport rate between two organs is
 * 
 * transport_rate = base_conductance * (total structural carbon)^m * (difference in mass fractions).
 * 
 * Important note: in his paper, Thornley uses a resistance rather than a conductance. We have chosen
 * to use a conductance instead to make it easy to close off a possible link between two organs. For
 * example, there is typically no direct link between the leaf and root of a plant. This link can be
 * closed by setting the `substrate_conductace_leaf_to_root` input to be zero. If we had used resistances
 * instead, we would have needed to specify an infinite resistance, which is not straightforward to do.
 */
class thornley_transport_calculator : public direct_module
{
   public:
    thornley_transport_calculator(
        std::vector<thornley_nutrient_dynamics::transport_link> const& organ_links,
        state_map const& input_quantities,
        state_map* output_quantities);

    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_name() { return "thornley_transport_calculator"; }

   private:
    std::vector<thornley_nutrient_dynamics::transport_link> const organ_links;

    // Pointers to input parameters
    std::vector<std::pair<double const*, double const*>> const substrate_carbon_ips;
    std::vector<std::pair<double const*, double const*>> const structural_carbon_ips;
    std::vector<double const*> const substrate_conductance_ips;

    // References to input parameters
    double const& transportation_beta_exponent;
    double const& Pod_start_dvi;
    double const& stop_growth_dvi;
    double const& DVI;

    // Pointers to output parameters
    std::vector<double*> const substrate_transport_ops;

   protected:
    void do_multi_organ_operation() const;
    static std::vector<std::string> get_inputs(std::vector<thornley_nutrient_dynamics::transport_link> const& organ_links);
    static std::vector<std::string> get_outputs(std::vector<thornley_nutrient_dynamics::transport_link> const& organ_links);
};

#endif
