# divEnv17032021.Rmd 
# libraries
library(glmulti)
# workspace
load(".RData")
# chunk023: species-by-species model for piAll
piAll_by_species_T4_4.glm <- list()
for (i in 1:7) 
{
  varList <- names(allDivEnvStar.df)[c(3:(4+nbLayers[i]))]
  piAll_by_species_T4_4.glm[[i]] <- 
    glmulti(y = "all_pi", xr = varList,
            data = bySpecies.ls[[i]], 
            level = 2, method = "h",fitfunction = lm, crit = 'bic', plotty = F)
}
names(piAll_by_species_T4_4.glm) <- names(bySpecies.ls)
# extra command to store output
save(piAll_by_species_T4_4.glm, file = "piAll_by_species_T4_4.glm.RData")
