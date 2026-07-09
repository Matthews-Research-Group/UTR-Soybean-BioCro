#include <vector>
#include <string>
#include "thornley_nutrient_dynamics_lsrp.h"

using thornley_nutrient_dynamics_lsrp::organ_list;
using thornley_nutrient_dynamics_lsrp::organ_transport_links;

////////////////////////////
// utilization calculator //
////////////////////////////

thornley_utilization_calculator_lsrp::thornley_utilization_calculator_lsrp(
    state_map const& input_quantities,
    state_map* output_quantities)
    : thornley_utilization_calculator(
          organ_list,
          input_quantities,
          output_quantities)
{
}

std::vector<std::string> thornley_utilization_calculator_lsrp::get_inputs()
{
    return thornley_utilization_calculator::get_inputs(organ_list);
}

std::vector<std::string> thornley_utilization_calculator_lsrp::get_outputs()
{
    return thornley_utilization_calculator::get_outputs(organ_list);
}

std::string thornley_utilization_calculator_lsrp::get_name()
{
    return "thornley_utilization_calculator_lsrp";
}

void thornley_utilization_calculator_lsrp::do_operation() const
{
    thornley_utilization_calculator::do_multi_organ_operation();
}

/////////////////
// utilization //
/////////////////

thornley_utilization_lsrp::thornley_utilization_lsrp(
    state_map const& input_quantities,
    state_map* output_quantities)
    : thornley_utilization(
          organ_list,
          input_quantities,
          output_quantities)
{
}

std::vector<std::string> thornley_utilization_lsrp::get_inputs()
{
    return thornley_utilization::get_inputs(organ_list);
}

std::vector<std::string> thornley_utilization_lsrp::get_outputs()
{
    return thornley_utilization::get_outputs(organ_list);
}

std::string thornley_utilization_lsrp::get_name()
{
    return "thornley_utilization_lsrp";
}

void thornley_utilization_lsrp::do_operation() const
{
    thornley_utilization::do_multi_organ_operation();
}

//////////////////////////
// transport calculator //
//////////////////////////

thornley_transport_calculator_lsrp::thornley_transport_calculator_lsrp(
    state_map const& input_quantities,
    state_map* output_quantities)
    : thornley_transport_calculator(
          organ_transport_links,
          input_quantities,
          output_quantities)
{
}

std::vector<std::string> thornley_transport_calculator_lsrp::get_inputs()
{
    return thornley_transport_calculator::get_inputs(organ_transport_links);
}

std::vector<std::string> thornley_transport_calculator_lsrp::get_outputs()
{
    return thornley_transport_calculator::get_outputs(organ_transport_links);
}

std::string thornley_transport_calculator_lsrp::get_name()
{
    return "thornley_transport_calculator_lsrp";
}

void thornley_transport_calculator_lsrp::do_operation() const
{
    thornley_transport_calculator::do_multi_organ_operation();
}
///////////////
// transport //
///////////////

thornley_transport_lsrp::thornley_transport_lsrp(
    state_map const& input_quantities,
    state_map* output_quantities)
    : thornley_transport(
          organ_transport_links,
          input_quantities,
          output_quantities)
{
}

std::vector<std::string> thornley_transport_lsrp::get_inputs()
{
    return thornley_transport::get_inputs(organ_transport_links);
}

std::vector<std::string> thornley_transport_lsrp::get_outputs()
{
    return thornley_transport::get_outputs(organ_transport_links);
}

std::string thornley_transport_lsrp::get_name()
{
    return "thornley_transport_lsrp";
}

void thornley_transport_lsrp::do_operation() const
{
    thornley_transport::do_multi_organ_operation();
}

///////////////////////
// biomass calculator//
///////////////////////

thornley_biomass_calculator_lsrp::thornley_biomass_calculator_lsrp(
    state_map const& input_quantities,
    state_map* output_quantities)
    : thornley_biomass_calculator(
          organ_list,
          input_quantities,
          output_quantities)
{
}

std::vector<std::string> thornley_biomass_calculator_lsrp::get_inputs()
{
    return thornley_biomass_calculator::get_inputs(organ_list);
}

std::vector<std::string> thornley_biomass_calculator_lsrp::get_outputs()
{
    return thornley_biomass_calculator::get_outputs(organ_list);
}

std::string thornley_biomass_calculator_lsrp::get_name()
{
    return "thornley_biomass_calculator_lsrp";
}


void thornley_biomass_calculator_lsrp::do_operation() const
{
    thornley_biomass_calculator::do_multi_organ_operation();
}