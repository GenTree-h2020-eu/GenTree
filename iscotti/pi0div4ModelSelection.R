# divEnv17032021.Rmd 
# libraries
library(glmulti)
# workspace
load(".RData")
# chunk024: species-by-species model for pi0div4
pi0div4_by_species_T4_4.glm <- list()
for (i in 1:7) 
{
  varList <- names(allDivEnvStar.df)[c(3:(4+nbLayers[i]))]
  piAll_by_species_T4_4.glm[[i]] <- 
    glmulti(y = "pi0div4", xr = varList,
            data = bySpecies.ls[[i]], 
            level = 2, method = "h",fitfunction = lm, crit = 'bic', plotty = F)
}
names(pi0div4_by_species_T4_4.glm) <- names(bySpecies.ls)
# extra command to store output
save(pi0div4_by_species_T4_4.glm, file = "pi0div4_by_species_T4_4.glm.RData")
