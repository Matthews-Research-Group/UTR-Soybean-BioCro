#ifndef THORNLEY_NUTRIENT_DYNAMICS_H
#define THORNLEY_NUTRIENT_DYNAMICS_H

#include <vector>
#include <string>
#include <unordered_map>

namespace thornley_nutrient_dynamics
{
/**
 * @class organ
 * @brief Allows the user to customize an organ by specifying an external source
 * of carbon substrate and/or designating it as a "storage organ," which may be
 * a source of carbon substrate rather than a sink. (Note: storage organ abilities
 * are not yet implemented.)
 */
class organ
{
   public:
    organ(
        std::string organ_name,
        std::string substrate_carbon_source,
        bool storage_organ)
        : organ_name(organ_name),
          substrate_carbon_source(substrate_carbon_source),
          storage_organ(storage_organ)
    {
    }
    std::string name() const { return organ_name; }
    bool has_external_carbon_substrate_source() const { return !substrate_carbon_source.empty(); }
    std::string external_substrate_carbon_source() const { return substrate_carbon_source; }
    bool is_storage_organ() const { return storage_organ; }

   private:
    std::string const organ_name;
    std::string const substrate_carbon_source;
    bool const storage_organ;
};

/** A library of available organs */
std::unordered_map<std::string, organ const> const organ_library = {
    {"FineRoot",        organ("FineRoot",       "",                         false)},
    {"Grain",           organ("Grain",          "",                         false)},
    {"Leaf",            organ("Leaf",           "canopy_assimilation_rate", false)},
    {"Rhizome",         organ("Rhizome",        "",                         true)},
    {"Root",            organ("Root",           "",                         false)},
    {"Seed",            organ("Seed",           "",                         true)},
    {"Stem",            organ("Stem",           "",                         false)},
    {"StructuralRoot",  organ("StructuralRoot", "",                         false)},
    {"Tuber",           organ("Tuber",          "",                         true)},
    {"Pod",             organ("Pod",            "",                         false)}// pod could have photosynthesis according to Zhang et al (2017)
                                                                                    // https://www.nature.com/articles/srep42026#:~:text=Most%20of%20the%20crop%20yield,is%20the%20result%20of%20photosynthesis.&text=Some%20non%2Dleaf%20organs%20contain,%2C9%2C10%2C11.
                                                                                    // "On a fresh weight basis, leaves photosynthesized 6.3–7.4 times more actively than pods. 
                                                                                    // However, on a chlorophyll basis, the pods had 1.5–1.8 times greater photosynthetic activity than the leaves." (Andrews & Svec, 1975)
};

/** A typedef for specifying transport links between organs */
using transport_link = std::pair<organ, organ>;

std::vector<std::string> generate_multi_organ_quantity_names(
    std::vector<organ> const& organs,
    std::vector<std::string> const& quantity_names);

std::vector<std::string> generate_quantity_names_from_transport_links(
    std::vector<transport_link> const& organ_links,
    std::vector<std::string> const& quantity_names);

std::vector<double const*> get_multi_organ_ips(
    state_map const& input_quantities,
    std::vector<organ> const& organs,
    std::string const& base_name);

std::vector<std::pair<double const*, double const*>> get_transport_link_ips(
    state_map const& input_quantities,
    std::vector<transport_link> const& organ_links,
    std::string const& base_name);

std::vector<double*> get_multi_organ_ops(
    state_map* output_quantities,
    std::vector<organ> const& organs,
    std::string const& base_name);

std::vector<std::pair<double*, double*>> get_transport_link_ops(
    state_map* output_quantities,
    std::vector<transport_link> const& organ_links,
    std::string const& base_name);

std::vector<std::string> get_external_substrate_quantity_names(
    std::vector<organ> const& organs);

std::vector<double const*> get_external_substrate_ips(
    state_map const& input_quantities,
    std::vector<organ> const& organs);

std::vector<std::string> generate_pairwise_names(std::string prefix, std::vector<transport_link> const& organ_links);

double hill_reaction_rate(
    double concentration,
    double hill_coefficient,
    double max_reaction_rate,
    double concentration_at_half_max);

// Define some constants used by the models
double const hill_coefficient = 2;

double senescence_logistic_fraction(
    double DVI,
    double max_senescence_fraction,
    double alpha,
    double beta);

}  // namespace thornley_nutrient_dynamics

#endif
