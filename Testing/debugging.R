# Preliminaries
library(BioCro)
library(PhotoGEA) # for pdf_print
library(lattice)

rm(list=ls())

# Choose settings
SAVE_TO_PDF <- TRUE

RUN_FULL_SIMULATIONS <- TRUE

# Define a soybean model for debugging
soybean_debug <- list(
    direct_modules = c(
        'BioCro:stomata_water_stress_linear',
        'BioCro:soybean_development_rate_calculator',
        'UTRSoybeanBML:thornley_utilization_calculator_lsrp',
        'UTRSoybeanBML:thornley_transport_calculator_lsrp',
        'UTRSoybeanBML:thornley_biomass_calculator_lsrp',
        'BioCro:parameter_calculator',
        'BioCro:soil_evaporation',
        'BioCro:solar_position_michalsky',
        'BioCro:shortwave_atmospheric_scattering',
        'BioCro:incident_shortwave_from_ground_par',
        'BioCro:ten_layer_canopy_properties',
        'BioCro:ten_layer_c3_canopy',
        'BioCro:ten_layer_canopy_integrator'
    ),
    differential_modules = c(
        'UTRSoybeanBML:thornley_utilization_lsrp',
        'UTRSoybeanBML:thornley_transport_lsrp',
        'UTRSoybeanBML:thornley_biomass_lsrp',
        'BioCro:two_layer_soil_profile',
        'BioCro:thermal_time_linear'
    ),
    ode_solver = list(
        type                               = 'boost_rkck54',
        output_step_size                   = 1.0,
        adaptive_rel_error_tol             = 1e-4,
        adaptive_abs_error_tol             = 1e-4,
        adaptive_max_steps                 = 200
    ),
    initial_values = list(
        Leaf_respiration_loss              = 0,
        Stem_respiration_loss              = 0,
        Root_respiration_loss              = 0,
        Pod_respiration_loss               = 0,
        Leaf_senescence_loss               = 0,
        Stem_senescence_loss               = 0,
        Root_senescence_loss               = 0,
        Pod_senescence_loss                = 0,
        TTc                                = 0,
        soil_water_content                 = 0.32,
        cws1                               = 0.32,
        cws2                               = 0.32,
        Rhizome                            = 1.00E-07,
        RhizomeLitter                      = 0,
        Leaf                               = 0.035505,
        Stem                               = 0.019725,
        Root                               = 0.02367,
        Pod                                = 0.001,
        Leaf_substrate_carbon              = 0.014202,
        Leaf_structural_carbon             = 0.127818,
        Stem_substrate_carbon              = 0.00789,
        Stem_structural_carbon             = 0.07101,
        Root_substrate_carbon              = 0.009468,
        Root_structural_carbon             = 0.085212,
        Pod_substrate_carbon               = 4.00E-04,
        Pod_structural_carbon              = 0.0036
    ),
    parameters = list(
        timestep                           = 1,
        Leaf_carbon_to_mass_factor         = 0.2442836,
        Leaf_utilization_rate_constant     = 0.02270617,
        Leaf_utilization_km                = 1.31624494,
        Leaf_respiration_factor            = 0,
        Stem_carbon_to_mass_factor         = 0.2442836,
        Stem_utilization_rate_constant     = 0.04000102,
        Stem_utilization_km                = 1.47867831,
        Stem_respiration_factor            = 0.08135216,
        Root_carbon_to_mass_factor         = 0.2442836,
        Root_utilization_rate_constant     = 0.03756,
        Root_utilization_km                = 1.4318081,
        Root_respiration_factor            = 0.08135216,
        Pod_carbon_to_mass_factor          = 0.2442836,
        Pod_utilization_rate_constant      = 0.55270866,
        Pod_utilization_km                 = 0.90805499,
        Pod_respiration_factor             = 0.08135216,
        substrate_conductance_Leaf_to_Stem = 0.685889,
        substrate_conductance_Stem_to_Root = 0.51311881,
        substrate_conductance_Stem_to_Pod  = 1.88708826,
        transportation_beta_exponent       = 1,
        Pod_start_dvi                      = 1,
        stop_growth_dvi                    = 2,
        Leaf_senescence_rate_max           = 0.06991942,
        Stem_senescence_rate_max           = 0.09900254,
        Root_senescence_rate_max           = 0.01206099,
        Pod_senescence_rate_max            = 0,
        Leaf_senescence_alpha              = 17.7650363,
        Stem_senescence_alpha              = 24.41046211,
        Root_senescence_alpha              = 6.97516144,
        Pod_senescence_alpha               = 110,
        Leaf_senescence_beta               = -9.15102399,
        Stem_senescence_beta               = -9.92872499,
        Root_senescence_beta               = -2.01417925,
        Pod_senescence_beta                = -0.4,
        Leaf_senescence_reuse_factor       = 0.6,
        Stem_senescence_reuse_factor       = 0.6,
        Root_senescence_reuse_factor       = 0.6,
        Pod_senescence_reuse_factor        = 0.6,
        iSp                                = 2.5,
        Sp_thermal_time_decay              = 0,
        LeafN_0                            = 2,
        LeafN                              = 2,
        vmax_n_intercept                   = 0,
        vmax1                              = 111.2,
        alphab1                            = 0,
        alpha1                             = 32.5,
        soil_air_entry                     = -2.6,
        soil_b_coefficient                 = 5.2,
        soil_bulk_density                  = 1.35,
        soil_clay_content                  = 0.34,
        soil_field_capacity                = 0.32,
        soil_sand_content                  = 0.32,
        soil_saturated_conductivity        = 6.40E-05,
        soil_saturation_capacity           = 0.52,
        soil_silt_content                  = 0.34,
        soil_wilting_point                 = 0.2,
        rsec                               = 0.2,
        soil_clod_size                     = 0.04,
        soil_reflectance                   = 0.2,
        soil_transmission                  = 0.01,
        specific_heat_of_air               = 1010,
        stefan_boltzman                    = 5.67E-08,
        maturity_group                     = 3,
        Tbase_emr                          = 10,
        TTemr_threshold                    = 60,
        Rmax_emrV0                         = 0.199,
        Tmin_emrV0                         = 5,
        Topt_emrV0                         = 31.5,
        Tmax_emrV0                         = 45,
        Tmin_R0R1                          = 5,
        Topt_R0R1                          = 31.5,
        Tmax_R0R1                          = 45,
        Tmin_R1R7                          = 0,
        Topt_R1R7                          = 21.5,
        Tmax_R1R7                          = 38.7,
        sowing_time                        = 0,
        par_energy_fraction                = 0.5,
        par_energy_content                 = 0.235,
        soil_depth1                        = 0,
        soil_depth2                        = 2.5,
        soil_depth3                        = 10,
        wsFun                              = 2,
        hydrDist                           = 0,
        rfl                                = 0.2,
        rsdf                               = 0.44,
        phi1                               = 0.01,
        phi2                               = 1.5,
        tbase                              = 10,
        lat                                = 40,
        longitude                          = -88,
        time_zone_offset                   = -6,
        atmospheric_pressure               = 101325,
        atmospheric_transmittance          = 0.85,
        atmospheric_scattering             = 0.3,
        absorptivity_par                   = 0.8,
        chil                               = 0.81,
        kd                                 = 0.7,
        heightf                            = 3,
        kpLN                               = 0,
        leaf_reflectance                   = 0.2,
        leaf_transmittance                 = 0.2,
        lnfun                              = 0,
        jmax                               = 195,
        jmax_mature                        = 195,
        sf_jmax                            = 0.2,
        electrons_per_carboxylation        = 4.5,
        electrons_per_oxygenation          = 5.25,
        tpu_rate_max                       = 13,
        Rd                                 = 1.28,
        Catm                               = 372,
        O2                                 = 210,
        b0                                 = 0.008,
        b1                                 = 10.6,
        Gs_min                             = 0.001,
        theta                              = 0.76,
        minimum_gbw                        = 0.08,
        windspeed_height                   = 5,
        beta_PSII                          = 0.5,
        growth_respiration_fraction        = 0
    )
)

# Load all the weather data
debugging_weather <- list(
    `2002` = read.csv('2002_drivers.csv'),
    `2004` = read.csv('2004_drivers.csv'),
    `2005` = read.csv('2005_drivers.csv'),
    `2006` = read.csv('2006_drivers.csv')
)

# Check the inputs for one year
with(soybean_debug, {validate_dynamical_system_inputs(
    initial_values,
    parameters,
    debugging_weather[['2002']],
    direct_modules,
    differential_modules,
    verbose = TRUE
)})

# Define helping function for running and plotting a year of results
run_year <- function(yn) {
    # Get the drivers
    drivers <- debugging_weather[[yn]]

    # Run the model
    biocro_result <- with(soybean_debug, {run_biocro(
        initial_values,
        parameters,
        drivers,
        direct_modules,
        differential_modules,
        ode_solver,
        verbose = TRUE
    )})

    # Save results
    write.csv(biocro_result, file = file.path('debug_outputs', paste0(yn, '_results.csv')), row.names = FALSE)

    # Find the first NA leaf value and time, if it exists
    na_indx <- which(is.na(biocro_result$Leaf))[1]
    na_time <- if (!is.na(na_indx)) {
        biocro_result$time[na_indx]
    } else {
        max(biocro_result$time)
    }

    # Get the time limits from the drivers
    drivers <- add_time_to_weather_data(drivers)

    # Find the first time DVI exceeds Pod_start_dvi
    pod_start_indx <- which(drivers$DVI > soybean_debug$parameters$Pod_start_dvi)[1]
    start_time <- if (!is.na(pod_start_indx)) {
        drivers$time[pod_start_indx]
    } else {
        min(drivers$time)
    }

    # Find the first time DVI exceeds stop_growth_dvi
    stop_time_indx <- which(drivers$DVI > soybean_debug$parameters$stop_growth_dvi)[1]
    stop_time <- if (!is.na(stop_time_indx)) {
        drivers$time[stop_time_indx]
    } else {
        min(drivers$time)
    }

    # Calculate "mass fractions"
    biocro_result <- within(biocro_result, {
        Leaf_mass_fraction = Leaf_substrate_carbon / Leaf
        Stem_mass_fraction = Stem_substrate_carbon / Leaf
        Root_mass_fraction = Root_substrate_carbon / Leaf
        Pod_mass_fraction = Pod_substrate_carbon / Leaf
    })

    # Calculate "proportional masses"
    biocro_result <- within(biocro_result, {
        Leaf_proportional_mass = Leaf_structural_carbon * soybean_debug$parameters$Leaf_carbon_to_mass_factor
        Stem_proportional_mass = Stem_structural_carbon * soybean_debug$parameters$Stem_carbon_to_mass_factor
        Root_proportional_mass = Root_structural_carbon * soybean_debug$parameters$Root_carbon_to_mass_factor
        Pod_proportional_mass = Pod_structural_carbon * soybean_debug$parameters$Pod_carbon_to_mass_factor
    })

    # Plot results
    time_lim <- c(min(drivers$time), max(drivers$time))
    time_lab <- paste0('Day of year (', yn, ')')

    pdf_print(
        xyplot(
            lai ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'LAI',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_lai.pdf'))
    )

    pdf_print(
        xyplot(
            Leaf + Stem + Root + Pod ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Biomass',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_biomass.pdf'))
    )

    pdf_print(
        xyplot(
            Leaf_mass_fraction + Stem_mass_fraction + Root_mass_fraction + Pod_mass_fraction ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Mass fraction (substrate_carbon / total mass)',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_mass_fraction.pdf'))
    )

    pdf_print(
        xyplot(
            Leaf_substrate_carbon + Stem_substrate_carbon + Root_substrate_carbon + Pod_substrate_carbon ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Substrate carbon',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_substrate_carbon.pdf'))
    )

    pdf_print(
        xyplot(
            Leaf_structural_carbon + Stem_structural_carbon + Root_structural_carbon + Pod_structural_carbon ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Structural carbon',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_structural_carbon.pdf'))
    )

    pdf_print(
        xyplot(
            Leaf_utilization_rate + Stem_utilization_rate + Root_utilization_rate + Pod_utilization_rate ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Utilization rate',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_utilization_rate.pdf'))
    )

    pdf_print(
        xyplot(
            Leaf_respiration_loss + Stem_respiration_loss + Root_respiration_loss + Pod_respiration_loss ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Respiration loss',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_respiration_loss.pdf'))
    )

    pdf_print(
        xyplot(
            Leaf_senescence_loss + Stem_senescence_loss + Root_senescence_loss + Pod_senescence_loss ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Senescence loss',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_senescence_loss.pdf'))
    )

    pdf_print(
        xyplot(
            Leaf_senescence_rate + Stem_senescence_rate + Root_senescence_rate + Pod_senescence_rate ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Senescence rate',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_senescence_rate.pdf'))
    )

    pdf_print(
        xyplot(
            Leaf_total_C_change_per_m2 + Stem_total_C_change_per_m2 + Root_total_C_change_per_m2 + Pod_total_C_change_per_m2 ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Total C change per m^2',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_total_C_change_per_m2.pdf'))
    )

    pdf_print(
        xyplot(
            substrate_transport_Leaf_to_Stem + substrate_transport_Stem_to_Pod + substrate_transport_Stem_to_Root ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Substrate transport',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_substrate_transport.pdf'))
    )

    DVI_step = 0.25

    pdf_print(
        xyplot(
            Leaf_total_C_change_per_m2  ~ Leaf,
            group = factor(DVI_step * floor(DVI / DVI_step)),
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            grid = TRUE,
            main = paste0(yn, '\nGrouped by DVI range'),
            par.settings = list(
                superpose.line = list(col = multi_curve_colors()),
                superpose.symbol = list(col = multi_curve_colors(), pch = 16)
            )
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_leaf_phase_space.pdf'))
    )

    pdf_print(
        xyplot(
            Leaf_proportional_mass + Leaf  ~ Leaf_structural_carbon,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            grid = TRUE,
            main = paste0(yn, '\nLeaf_proportional_mass = Leaf_structural_carbon * Leaf_carbon_to_mass_factor'),
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_leaf_mass_space.pdf'))
    )

    pdf_print(
        xyplot(
            TTc ~ time,
            data = biocro_result,
            type = 'l',
            auto = TRUE,
            xlim = time_lim,
            xlab = time_lab,
            ylab = 'Thermal time',
            main = yn,
            panel = function(...) {
                panel.xyplot(...)
                panel.lines(c(100, -100) ~ c(na_time, na_time), lty = 2, col = 'darkgray')
                panel.lines(c(100, -100) ~ c(stop_time, stop_time), lty = 1, col = 'black')
                panel.lines(c(100, -100) ~ c(start_time, start_time), lty = 4, col = 'red')
            }
        ),
        width = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path('debug_outputs', paste0(yn, '_thermal_time.pdf'))
    )

    return(biocro_result)
}

if (RUN_FULL_SIMULATIONS) {
    results <- lapply(names(debugging_weather), run_year)
}
