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


mod.iid.somm <- mmer(fixed = y ~ factweek + variety + block + factweek:variety + factweek:block, 
                     rcov  = ~ units, # default = iid
                     data  = dat, verbose=F)

# Extract variance component estimates
mod.iid.somm.VC <- summary(mod.iid.somm)$varcomp 




sasvc <- readr::read_delim("SAS/autoregressive_time_series_results_VC.txt", delim="\t")

sasvciid <- sasvc %>%
  filter(mod=="iid")



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


fixed.rho <- 0.7 # random guess

mod.ar1.somm <- mmer(fixed  = y ~ factweek + variety + block + factweek:variety + factweek:block, 
                     random = ~ vs(factweek, Gu = AR1(factweek, rho = fixed.rho)),
                     rcov   = ~ units, # default = iid
                     data   = dat, verbose=F)

# Extract variance component estimates
mod.ar1.somm.VC <- summary(mod.ar1.somm)$varcomp 




sasvcar1 <- sasvc %>%
  filter(mod=="ar1")



AIC.nlme <- aictab(list(mod.iid.nlme, mod.ar1.nlme), modnames = c("iid", "ar1"))


AIC.glmm <- aictab(list(mod.iid.glmm, mod.ar1.glmm), modnames = c("iid", "ar1"))


somm.mods <- list(mod.iid.somm, mod.ar1.somm)

AIC.somm <- tibble(
  Modnames = c("iid", "ar1"),
  AIC = somm.mods %>% map("AIC") %>% unlist,
  LL  = somm.mods %>% map("monitor") %>% lapply(. %>% `[`(1, ncol(.))) %>% unlist) %>% # last element in first row of "monitor"
  arrange(AIC)


sasaic <- readr::read_delim("SAS/autoregressive_time_series_results_AIC.txt", delim="\t") %>% 
  mutate(Descr = str_remove(Descr, "\\s*\\([^\\)]+\\)") %>% str_trim) %>% 
  pivot_wider(names_from=Descr, values_from=Value)



tibble(mod="iid", estimate="variance") %>%
  cbind(., 
        mod.iid.nlme.VC %>% dplyr::select(Variance) %>% rename(nlme=Variance),
        mod.iid.nlme.VC %>% mutate(lme4=NA) %>% dplyr::select(lme4),
        mod.iid.glmm.VC %>% filter(effect=="Intercept") %>% dplyr::select(variance) %>% rename(glmmTMB=variance),
        mod.iid.somm.VC %>% dplyr::select(VarComp) %>% rename(sommer=VarComp),
        sasvc %>% filter(mod=="iid") %>% dplyr::select(Estimate) %>% rename(SAS=Estimate)) %>% 
  magrittr::set_rownames(NULL) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

tibble(mod=c("ar1", "ar1"), estimate=c("variance", "correlation")) %>%
  cbind(., 
        mod.ar1.nlme.VC %>% dplyr::select(Variance, rho) %>% 
          pivot_longer(names_to = "estimate", values_to = "nlme", cols = 1:2) %>% 
          dplyr::select(nlme),
        tibble(lme4=c(NA, NA)),
        tibble(glmmTMB=
                 c(mod.ar1.glmm.VC %>%
                     pluck("Variance Components") %>% 
                     filter(effect=="factweek1") %>% 
                     dplyr::select(variance), 
                   mod.ar1.glmm.VC %>%
                     pluck("Cor") %>% pluck("plot") %>%  
                     pluck(2,1))%>% unlist),
        tibble(sommer=c(NA, NA)),
        sasvc %>% 
          arrange(desc(CovParm)) %>% 
          filter(mod=="ar1") %>% 
          dplyr::select(Estimate) %>% 
          rename(SAS=Estimate)) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

tibble(mod=c("iid", "ar1"), estimate=c("AIC", "AIC")) %>%
  cbind(.,
        AIC.nlme %>% dplyr::select(AICc) %>% rename(nlme=AICc),
        tibble(lme4=c(NA, NA)),
        AIC.glmm %>% dplyr::select(AICc) %>% rename(glmmTMB=AICc),
        tibble(sommer=c(AIC.somm %>% filter(Modnames=="iid") %>% pull(AIC), NA)),
        sasaic %>% dplyr::select(AIC) %>% rename(SAS=AIC)) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)
