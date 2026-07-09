#ifndef THORNLEY_NUTRIENT_DYNAMICS_LSRP_H
#define THORNLEY_NUTRIENT_DYNAMICS_LSRP_H

#include <vector>
#include <string>
#include "../framework/state_map.h"
#include "thornley_nutrient_dynamics.h"
#include "thornley_utilization_calculator.h"
#include "thornley_utilization.h"
#include "thornley_transport_calculator.h"
#include "thornley_transport.h"
#include "thornley_biomass_calculator.h"

using thornley_nutrient_dynamics::organ;
using thornley_nutrient_dynamics::organ_library;
using thornley_nutrient_dynamics::transport_link;

namespace thornley_nutrient_dynamics_lsrp
{
/** Specify the organs to use for the `lsrg` version */
std::vector<organ> const organ_list = {
    organ_library.at("Leaf"),  // `l` is for Leaf
    organ_library.at("Stem"),  // `s` is for Stem
    organ_library.at("Root"),  // `r` is for Root
    organ_library.at("Pod")    // `p` is for Pod
};

/** Specify the substrate transport links to use for the `lsrg` version */
std::vector<transport_link> const organ_transport_links = {
    transport_link(organ_library.at("Leaf"), organ_library.at("Stem")),  // Leaf to Stem
    transport_link(organ_library.at("Stem"), organ_library.at("Root")),  // Stem to Root
    transport_link(organ_library.at("Stem"), organ_library.at("Pod"))    // Stem to Pod
};

}  // namespace thornley_nutrient_dynamics_lsrp

/**
 * @class thornley_utilization_calculator_lsrp
 * 
 * @brief A child class of thornley_utilization_calculator where the organs have been set to
 * `leaf`, `stem`, `root`, and `pod` (hence the name `lsrp`). Instances of this class can be created using the
 * module wrapper factory.
 */
class thornley_utilization_calculator_lsrp : public thornley_utilization_calculator
{
   public:
    thornley_utilization_calculator_lsrp(
        state_map const& input_quantities,
        state_map* output_quantities);

    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_name();


   private:
    // Main operation
    void do_operation() const;
};

/**
 * @class thornley_utilization_lsrp
 * 
 * @brief A child class of thornley_utilization where the organs have been set to
 * `leaf`, `stem`, `root`, and `pod` (hence the name `lsrp`). Instances of this class can be created using the
 * module wrapper factory.
 */
class thornley_utilization_lsrp : public thornley_utilization
{
   public:
    thornley_utilization_lsrp(
        state_map const& input_quantities,
        state_map* output_quantities);

    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_name();

   private:
    // Main operation
    void do_operation() const;
};

/**
 * @class thornley_transport_calculator_lsrp
 * 
 * @brief A child class of thornley_transport_calculator where the organs have been set to
 * `leaf`, `stem`, `root`, and `pod` (hence the name `lsrp`). Instances of this class can be created using the
 * module wrapper factory.
 */
class thornley_transport_calculator_lsrp : public thornley_transport_calculator
{
   public:
    thornley_transport_calculator_lsrp(
        state_map const& input_quantities,
        state_map* output_quantities);

    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_name();

   private:
    // Main operation
    void do_operation() const;
};

/**
 * @class thornley_transport_lsrp
 * 
 * @brief A child class of thornley_transport where the organs have been set to
 * `leaf`, `stem`, `root`, and `pod` (hence the name `lsrp`). Instances of this class can be created using the
 * module wrapper factory.
 */
class thornley_transport_lsrp : public thornley_transport
{
   public:
    thornley_transport_lsrp(
        state_map const& input_quantities,
        state_map* output_quantities);

    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_name();

   private:
    // Main operation
    void do_operation() const;
};
/**
 * @class thornley_biomass_calculator_lsrp
 * 
 * @brief A child class of thornley_biomass_calculator where the organs have been set to
 * `leaf`, `stem`, `root`, and `pod` (hence the name `lsrp`). Instances of this class can be created using the
 * module wrapper factory.
 */

class thornley_biomass_calculator_lsrp : public thornley_biomass_calculator
{
   public:
    thornley_biomass_calculator_lsrp(
        state_map const& input_quantities,
        state_map* output_quantities);

    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_name();

   private:
    // Main operation
    void do_operation() const;
};

#endif