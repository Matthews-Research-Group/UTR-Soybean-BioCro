#include <vector>
#include <set>
#include <string>
#include <cmath>                                     // for pow and abs
#include <algorithm>                                 // for max
#include "../framework/module_helper_functions.h"    // for add_class_prefix_to_quantity_name, get_ip, get_op
#include "../framework/validate_dynamical_system.h"  // for string_set_to_string_vector
#include "thornley_nutrient_dynamics.h"

using thornley_nutrient_dynamics::organ;
using thornley_nutrient_dynamics::transport_link;

/**
 * @brief Adds organ name prefixes to the quantity names indicating that each
 * separate version corresponds to a different organ.
 */
std::vector<std::string> thornley_nutrient_dynamics::generate_multi_organ_quantity_names(
    std::vector<organ> const& organs,
    std::vector<std::string> const& quantity_names)
{
    std::vector<std::string> multi_organ_vector;
    for (organ const& x : organs) {
        for (std::string const& n : quantity_names) {
            if (n == ""){
                multi_organ_vector.push_back(x.name());
            }else{
                multi_organ_vector.push_back(add_class_prefix_to_quantity_name(x.name(), n));
            }
        }
    }
    return multi_organ_vector;
}

/**
 * @brief Adds organ name prefixes to the quantity names indicating that each
 * separate version corresponds to one of the organs included in the transport
 * links.
 */
std::vector<std::string> thornley_nutrient_dynamics::generate_quantity_names_from_transport_links(
    std::vector<transport_link> const& organ_links,
    std::vector<std::string> const& quantity_names)
{
    std::set<std::string> all_quantity_names;  // use a set to prevent duplicates
    for (transport_link const& x : organ_links) {
        for (std::string const& n : quantity_names) {
            if (n == ""){
                all_quantity_names.insert(x.first.name());
                all_quantity_names.insert(x.second.name());
            }else{
                all_quantity_names.insert(add_class_prefix_to_quantity_name(x.first.name(), n));
                all_quantity_names.insert(add_class_prefix_to_quantity_name(x.second.name(), n));
            }
        }
    }
    return string_set_to_string_vector(all_quantity_names);
}

/**
 * @brief Returns a vector of `input pointers` using quantity names determined by
 * the base name and the set of organs.
 */
std::vector<double const*> thornley_nutrient_dynamics::get_multi_organ_ips(
    state_map const& input_quantities,
    std::vector<organ> const& organs,
    std::string const& base_name)
{
    std::vector<std::string> quantity_names = thornley_nutrient_dynamics::generate_multi_organ_quantity_names(organs, {base_name});
    std::vector<double const*> multi_organ_references(quantity_names.size());
    for (size_t i = 0; i < quantity_names.size(); i++) {
        multi_organ_references[i] = get_ip(input_quantities, quantity_names[i]);
    }
    return multi_organ_references;
}

/**
 * @brief Returns a vector of pairs of `input pointers` using quantity names
 * determined by the base name and the pair of organs involved in the transport
 * link.
 */
std::vector<std::pair<double const*, double const*>> thornley_nutrient_dynamics::get_transport_link_ips(
    state_map const& input_quantities,
    std::vector<transport_link> const& organ_links,
    std::string const& base_name)
{
    std::vector<std::pair<double const*, double const*>> pointer_pairs;
    for (size_t i = 0; i < organ_links.size(); ++i) {
        std::string first_name;
        std::string second_name;
        if (base_name == ""){
            first_name = organ_links[i].first.name();
            second_name = organ_links[i].second.name();
        }else{
            first_name = add_class_prefix_to_quantity_name(organ_links[i].first.name(), base_name);
            second_name = add_class_prefix_to_quantity_name(organ_links[i].second.name(), base_name);
        }
        
        pointer_pairs.push_back(std::pair<double const*, double const*>(
            get_ip(input_quantities, first_name),
            get_ip(input_quantities, second_name)));
    }
    return pointer_pairs;
}

/**
 * @brief Returns a vector of `output pointers` using quantity names determined by
 * the base name and the set of organs.
 */
std::vector<double*> thornley_nutrient_dynamics::get_multi_organ_ops(
    state_map* output_quantities,
    std::vector<organ> const& organs,
    std::string const& base_name)
{
    std::vector<std::string> quantity_names = thornley_nutrient_dynamics::generate_multi_organ_quantity_names(organs, {base_name});
    std::vector<double*> multi_organ_references(quantity_names.size());
    for (size_t i = 0; i < quantity_names.size(); i++) {
        multi_organ_references[i] = get_op(output_quantities, quantity_names[i]);
    }
    return multi_organ_references;
}

/**
 * @brief Returns a vector of pairs of `output pointers` using quantity names
 * determined by the base name and the pair of organs involved in the transport
 * link.
 */
std::vector<std::pair<double*, double*>> thornley_nutrient_dynamics::get_transport_link_ops(
    state_map* output_quantities,
    std::vector<transport_link> const& organ_links,
    std::string const& base_name)
{
    std::vector<std::pair<double*, double*>> pointer_pairs;
    for (size_t i = 0; i < organ_links.size(); ++i) {
        std::string first_name = add_class_prefix_to_quantity_name(organ_links[i].first.name(), base_name);
        std::string second_name = add_class_prefix_to_quantity_name(organ_links[i].second.name(), base_name);
        pointer_pairs.push_back(std::pair<double*, double*>(
            get_op(output_quantities, first_name),
            get_op(output_quantities, second_name)));
    }
    return pointer_pairs;
}

/**
 * @brief Gets the names of any external carbon substrate sources required by the organs.
 */
std::vector<std::string> thornley_nutrient_dynamics::get_external_substrate_quantity_names(
    std::vector<organ> const& organs)
{
    std::vector<std::string> quantity_names;
    for (organ const& x : organs) {
        if (x.has_external_carbon_substrate_source()) {
            quantity_names.push_back(x.external_substrate_carbon_source());
        }
    }
    return quantity_names;
}

/**
 * @brief Gets a pointer to each organ's external carbon substrate source. If an organ
 * lacks an exernal source of carbon substrate, this function will return a NULL pointer
 * for the element corresponding to that organ.
 */
std::vector<double const*> thornley_nutrient_dynamics::get_external_substrate_ips(
    state_map const& input_quantities,
    std::vector<organ> const& organs)
{
    std::vector<double const*> ip_vector(organs.size());
    for (size_t i = 0; i < organs.size(); ++i) {
        if (organs[i].has_external_carbon_substrate_source()) {
            ip_vector[i] = get_ip(input_quantities, organs[i].external_substrate_carbon_source());
        } else {
            ip_vector[i] = NULL;
        }
    }
    return ip_vector;
}

/**
 * @brief A function that generates quantity names for quantities involved in the transport of substrate between pairs of organs.
 */
std::vector<std::string> thornley_nutrient_dynamics::generate_pairwise_names(std::string prefix, std::vector<transport_link> const& organ_links)
{
    std::vector<std::string> transport_names;
    for (size_t i = 0; i < organ_links.size(); ++i) {
        std::string quantity_name = prefix + "_" + organ_links[i].first.name() + "_to_" + organ_links[i].second.name();
        transport_names.push_back(quantity_name);
    }
    return transport_names;
}

/**
 * @brief A function that determines a reaction rate using the Hill equation, where the output
 * has the same units as the `max_reaction_rate` input. Note that when `hill_coefficient` is 1,
 * the Hill equation reduces to the Michaelis-Menten equation. See the Wikipedia page for more
 * details: https://en.wikipedia.org/wiki/Hill_equation_(biochemistry).
 * If the substrate concentration is negative, then 
 */
double thornley_nutrient_dynamics::hill_reaction_rate(
    double concentration,
    double hill_coefficient,
    double max_reaction_rate,
    double concentration_at_half_max)
{
    return max_reaction_rate / (1.0 + pow(concentration_at_half_max / std::max(1e-10, concentration), hill_coefficient));
}
/**
 * @brief A function that determines a reaction rate using the logistic equation, where the output
 * has the same units as the `max_reaction_rate` input. The logistic equations follow the formulation
 * that in Soybean-BioCro. 
 * Noted: The DVI is replaced with time - 210 for now since DVI cannot work properly yet.
 * Reference: Megan L Matthews, Amy Marshall-Colón, Justin M McGrath, Edward B Lochocki, Stephen P Long, 
 * Soybean-BioCro: a semi-mechanistic model of soybean growth, in silico Plants, Volume 4, Issue 1, 
 * 2022, diab032, https://doi.org/10.1093/insilicoplants/diab032
 */
double thornley_nutrient_dynamics::senescence_logistic_fraction(
    double DVI,
    double max_senescence_fraction,
    double alpha,
    double beta)
{   
    return  max_senescence_fraction / (1.0 + exp(alpha * (beta - DVI)));
}