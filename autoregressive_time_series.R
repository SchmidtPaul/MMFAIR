params <-
list(hilang = "sas")

# packages
pacman::p_load(dplyr, purrr, tibble, tidyr, stringr, # data handling
               ggplot2, viridis,            # plot
               nlme, lme4, glmmTMB, sommer, # mixed modelling
               AICcmodavg, mixedup)         # mixed model extractions

# data
dat <- agriTutorial::sorghum %>%
  rename(block = Replicate, plot = factplot) %>% 
  dplyr::select(y, variety, block, plot, factweek, varweek) %>% 
  as_tibble()




dat.wk1 <- dat %>% filter(factweek == "1") # subset data from first week only

mod.wk1 <- lm(formula = y ~ variety + block,
              data = dat.wk1)

mod.iid.nlme <- nlme::gls(model = y ~ factweek * (variety + block),
                          correlation = NULL, # default, i.e. homoscedastic, independent errors
                          data = dat)

# Extract variance component estimates
mod.iid.nlme.VC <- tibble(varstruct = "iid") %>% 
  mutate(sigma    = mod.iid.nlme$sigma) %>% 
  mutate(Variance = sigma^2)


dat <- dat %>%
  mutate(unit = 1:n() %>% as.factor) # new column with running number

mod.iid.glmm <- glmmTMB(formula = y ~ factweek * (variety + block) 
                        + (1 | unit),      # add random unit term to mimic error variance
                        dispformula = ~ 0, # fix original error variance to 0
                        REML = TRUE,       # needs to be stated since default = ML
                        data = dat) 

# Extract variance component estimates
mod.iid.glmm.VC <- mod.iid.glmm %>%
  extract_vc(ci_scale = "var")




mod.ar1.nlme <- nlme::gls(model = y ~ factweek * (variety + block),
                          correlation = corAR1(form = ~ varweek | plot),
                          data = dat)

# Extract variance component estimates
mod.ar1.nlme.VC <- tibble(varstruct = "ar(1)") %>%
  mutate(sigma    = mod.ar1.nlme$sigma,
         rho      = coef(mod.ar1.nlme$modelStruct$corStruct, unconstrained = FALSE)) %>%
  mutate(Variance = sigma^2,
         Corr1wk  = rho,
         Corr2wks = rho^2,
         Corr3wks = rho^3,
         Corr4wks = rho^4)

mod.ar1.glmm <- glmmTMB(formula = y ~ factweek * (variety + block)
                        + ar1(factweek + 0 | plot), # add ar1 structure as random term to mimic error variance
                        dispformula = ~ 0, # fix original error variance to 0
                        REML = TRUE,       # needs to be stated since default = ML
                        data = dat) 

# Extract variance component estimates
mod.ar1.glmm.VC <- mod.ar1.glmm %>%
  extract_vc(ci_scale = "var", show_cor = TRUE)


mod.iid.somm <- mmer(fixed  = y ~ factweek + variety + block + factweek:variety + factweek:block, 
                     rcov   = ~ vs(AR1(plot),factweek), 
                     data   = dat, verbose=F)
