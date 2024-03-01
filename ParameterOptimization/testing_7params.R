library(BioCro)
library(DEoptim)
library(lattice)
library(lhs)
library(ppcor)
# Clear workspace
rm(list=ls())

# Set work directory
setwd("C://Users//pxm72//Documents//Projects//Soybean-BioCro//ParameterOptimization")

source('multiyear_partitioning_senescence_coefs_optimization.R')

test_solver <- match.fun(soybean_optsolver[[1]])

# parameter upper limit
upperlim<-c(100,5e-3, 100, 0.1, 
            5e-3, 5e-2, 0.2)

# parameter lower limit
lowerlim<-c(10, 1e-5, 10, 0, 
            1e-3, 1e-5, 0.01)

difflim = upperlim - lowerlim

sample_size = 400

numofparams = length(upperlim)

lhs_matrix_7 = randomLHS(sample_size,numofparams)

for (i in 1:numofparams){
  lhs_matrix_7[,i] = lhs_matrix_7[,i]*difflim[i]+lowerlim[i]
}
lhs_matrix <- list()
lhs_matrix = matrix(numeric(sample_size*16) ,nrow = sample_size, ncol = 16)
lhs_matrix = data.frame(lhs_matrix)
colnames(lhs_matrix) = append(arg_names, "total_biomass")
# carbon to mass factor
lhs_matrix[,c(1,5,9)] = lhs_matrix_7[,1]

# utilization rate constant
lhs_matrix[,2] = 0.4* lhs_matrix_7[,2]
lhs_matrix[,6] = lhs_matrix_7[,2] # picked
lhs_matrix[,10] = 0.1* lhs_matrix_7[,2]

# utilization km
lhs_matrix[,3] = lhs_matrix_7[,3] # picked
lhs_matrix[,7] = 1.2 * lhs_matrix_7[,3]
lhs_matrix[,11] = 0.5* lhs_matrix_7[,3]

# respiration factor
lhs_matrix[,4] = 0 # Leaf respiration is accounted for by the canopy photosynthesis module
lhs_matrix[,8] = lhs_matrix_7[,4] # picked
lhs_matrix[,12] = 1.5* lhs_matrix_7[,4]

# substrate conductance and utilization beta component
lhs_matrix[,c(13,14,15)] = lhs_matrix_7[,c(5,6,7)]

lhs_matrix[,'total_biomass'] = 0

for (x in 1:sample_size){
  # initialize parameters being fitted as random values from a uniform distribution
  # soybean_solver <- partial_run_biocro(initial_state,
  #                                      parameters,
  #                                      weather.growingseason,
  #                                      steady_state_module_names,
  #                                      derivative_module_names,
  #                                      solver,
  #                                      arg_names,
  #                                      verbose = FALSE)
  # test_solver = match.fun(soybean_solver)
  result <- test_solver(lhs_matrix[x,1:15])
  last_row = tail(result, n = 1)
  lhs_matrix[x,'total_biomass'] = last_row$Stem + last_row$Leaf + last_row$Root
  # gc()
}
png(file=paste("lhs_7param_" , sample_size, "_72_prcc_pv.png", sep=""), width=1000, height=800)
attach(lhs_matrix)
par(mfrow=c(4,4))
par(mar=c(2,2,2,2))
plot(Leaf_structural_carbon_to_mass_factor, total_biomass, main = "Leaf structural carbon to mass factor")
plot(Leaf_utilization_rate_constant, total_biomass, main = "Leaf utilization rate constant")
plot(Leaf_utilization_km, total_biomass, main = "Leaf utilization km")
plot(Leaf_respiration_factor, total_biomass, main = "Leaf respiration factor")
plot(Stem_structural_carbon_to_mass_factor, total_biomass, main = "Stem structural carbon to mass factor")
plot(Stem_utilization_rate_constant,total_biomass, main = "Stem utilization rate constant")
plot(Stem_utilization_km, total_biomass, main = "Stem utilization km")
plot(Stem_respiration_factor, total_biomass, main = "Stem respiration factor")
plot(Root_structural_carbon_to_mass_factor, total_biomass, main = "Root structural carbon to mass factor")
plot(Root_utilization_rate_constant, total_biomass, main = "Root_ utilization rate constant")
plot(Root_utilization_km, total_biomass, main = "Root utilization km")
plot(Root_respiration_factor, total_biomass, main = "Root_respiration_factor")
plot(substrate_conductance_Leaf_to_Stem, total_biomass, main =  "substrate conductance Leaf to Stem")
plot(substrate_conductance_Stem_to_Root, total_biomass, main = "Substrate conductance Stem to Root")
plot(utilization_beta_exponent, total_biomass, main = "utilization beta exponent")
dev.off()
#lhs_matrix_no0 = lhs_matrix[which(lhs_matrix$total_biomass<1),]
# thornley.mlr = lm(total_biomass ~ Leaf_structural_carbon_to_mass_factor+
#                        Leaf_utilization_rate_constant+
#                        Leaf_utilization_km+
#                        Leaf_respiration_factor+
#                        Stem_structural_carbon_to_mass_factor+
#                        Stem_utilization_rate_constant+
#                        Stem_utilization_km+
#                        Stem_respiration_factor+
#                        Root_structural_carbon_to_mass_factor+
#                        Root_utilization_rate_constant+
#                        Root_utilization_km+
#                        Root_respiration_factor+
#                        substrate_conductance_Leaf_to_Stem+
#                        substrate_conductance_Stem_to_Root+
#                        utilization_beta_exponent, data=lhs_matrix)
#summary(thornley.mlr)

pcor_result = pcor(lhs_matrix, method = "spearman")
pcor_significant_pv = pcor_result$p.value[16, pcor_result$p.value[-16,16]<0.05]
lhs_matrix_prcc_pv = rbind(lhs_matrix, pcor_result$estimate[16,], pcor_result$p.value[16,])
save(lhs_matrix_prcc_pv,file = paste("lhs_7param_" , sample_size, "_72_prcc_pv.Rdata", sep=""))

# lhs_matrix_bi = lhs_matrix
# lhs_matrix_bi[,'lived'] = 1
# lhs_matrix_bi[lhs_matrix$total_biomass<1e-1,'lived'] = 0
# lhs_matrix_bi
# lhs_matrix_bi.mlr = lm(lived ~ Leaf_structural_carbon_to_mass_factor+
#                                               Leaf_utilization_rate_constant+
#                                               Leaf_utilization_km+
#                                               Leaf_respiration_factor+
#                                               Stem_structural_carbon_to_mass_factor+
#                                               Stem_utilization_rate_constant+
#                                               Stem_utilization_km+
#                                               Stem_respiration_factor+
#                                               Root_structural_carbon_to_mass_factor+
#                                               Root_utilization_rate_constant+
#                                               Root_utilization_km+
#                                               Root_respiration_factor+
#                                               substrate_conductance_Leaf_to_Stem+
#                                               substrate_conductance_Stem_to_Root+
#                                               utilization_beta_exponent, data=lhs_matrix_bi)
# 
# summary(lhs_matrix_bi.mlr)