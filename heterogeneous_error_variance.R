params <-
list(hilang = "sas")

options(knitr.kable.NA = ' ')
# packages for better formatting tables for html output
library(kableExtra)
library(formattable)
library(htmltools)
library(stringr)

# packages
pacman::p_load(dplyr, purrr, tibble, tidyr, # data handling
               nlme, lme4, glmmTMB, sommer, # mixed modelling
               AICcmodavg, broom.mixed)     # mixed model extractions
# data
dat <- agridat::mcconway.turnip %>%
  as_tibble() %>%
  mutate(densf = density %>% as.factor)
dat %>%
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE) %>% 
  scroll_box(height = "200px")

boxplot(yield ~ date, data = dat)

boxplot(yield ~ density, data = dat)



dat <- dat %>%
  mutate(date_densf = interaction(date, densf)) # needed for mod5

mod1.nlme <- nlme::lme(fixed   = yield ~ gen * date * densf,  
                       random  = ~ 1 | block, 
                       weights = NULL, # default, i.e. homoscedastic errors
                       data    = dat)

mod2.nlme <- mod1.nlme %>%
  update(weights = varIdent(form =  ~ 1 | date))

mod3.nlme <- mod1.nlme %>%
  update(weights = varIdent(form =  ~ 1 | densf))

mod4.nlme <- mod1.nlme %>%
  update(weights = varComb(varIdent(form =  ~ 1 | date),
                           varIdent(form =  ~ 1 | densf))) 

mod5.nlme <- mod1.nlme %>%
  update(weights = varIdent(form =  ~ 1 | date_densf))

mod1.lme4 <- lmer(formula = yield ~ gen * date * densf + (1 | block),
                  data    = dat)

# mod2.lme4 - not possible

# mod3.lme4 - not possible

# mod4.lme4 - not possible

# mod5.lme4 - not possible

dat <- dat %>%
  mutate(unit = 1:n() %>% as.factor) # new column with running number

mod1.glmm <- glmmTMB(formula = yield ~ gen * date * densf + (1 | block), 
                     dispformula = ~ 1,  # = default i.e. homoscedastic error variance
                     REML = TRUE,        # needs to be stated since default = ML
                     data = dat) 

mod1b.glmm <- mod1.glmm %>%
  update(. ~ . + (1 | unit), # add random term to mimic homoscedastic error variance
         dispformula =  ~ 0) # fix original error variance to 0

mod2.glmm <- mod1.glmm %>% 
  update(.~. + diag(date+0|unit), dispformula=~0)

mod3.glmm <- mod1.glmm %>% 
  update(.~. + diag(densf+0|unit), dispformula=~0)

# mod4/multiplicative variance structure possible? in progress.

mod5.glmm <- mod1.glmm %>% 
  update(.~. + diag(date:densf+0|unit), dispformula=~0)

mod1.glmm %>% VarCorr()

mod1b.glmm %>% VarCorr()

mod1.somm <- mmer(fixed  = yield ~ gen * date * densf, 
                  random = ~ block,
                  rcov   = ~ units, # default
                  data   = dat, verbose=F)

mod2.somm <- mmer(fixed  = yield ~ gen * date * densf, 
                  random = ~ block,
                  rcov   = ~ vs(ds(date), units), 
                  data   = dat, verbose=F)

mod3.somm <- mmer(fixed  = yield ~ gen * date * densf, 
                  random = ~ block,
                  rcov   = ~ vs(ds(densf), units), 
                  data   = dat, verbose=F)

# mod4/multiplicative variance structure possible? in progress.

mod5.somm <- mmer(fixed  = yield ~ gen * date * densf, 
                  random = ~ block,
                  rcov   = ~ vs(ds(date_densf), units), 
                  data   = dat, verbose=F)
summary(mod5.somm)$varcomp

mod5.somm <- mmer(fixed  = yield ~ gen * date * densf, 
                   random = ~ block,
                   rcov   = ~ vs(ds(date), ds(densf), units), 
                   data   = dat, verbose=F)
summary(mod5b.somm)$varcomp

mod1.nlme.VC <- tibble(grp = "homoscedastic", varStruct = 1) %>% 
  mutate(sigma         = mod1.nlme$sigma) %>% 
  mutate(StandardError = sigma*varStruct) %>% 
  mutate(Variance      = StandardError^2)


mod2.nlme.VC <- mod2.nlme$modelStruct$varStruct %>%
  coef(unconstrained = FALSE, allCoef = TRUE) %>%
  enframe(name = "grp", value = "varStruct") %>%
  mutate(sigma         = mod2.nlme$sigma) %>%
  mutate(StandardError = sigma * varStruct) %>%
  mutate(Variance      = StandardError ^ 2)


mod3.nlme.VC <- mod3.nlme$modelStruct$varStruct %>%
  coef(unconstrained = FALSE, allCoef = TRUE) %>%
  enframe(name = "grp", value = "varStruct") %>%
  mutate(sigma         = mod3.nlme$sigma) %>%
  mutate(StandardError = sigma * varStruct) %>%
  mutate(Variance      = StandardError ^ 2)


mod4.nlme.varStruct.A <- mod4.nlme$modelStruct$varStruct$A %>%
  coef(unconstrained = FALSE, allCoef = TRUE) %>%
  enframe(name = "grpA", value = "varStructA")

mod4.nlme.varStruct.B <- mod4.nlme$modelStruct$varStruct$B %>%
  coef(unconstrained = FALSE, allCoef = TRUE) %>%
  enframe(name = "grpB", value = "varStructB")

mod4.nlme.VC <- expand.grid(mod4.nlme.varStruct.A$grpA,
                            mod4.nlme.varStruct.B$grpB,
                            stringsAsFactors = FALSE) %>%
  rename(grpA = Var1, grpB = Var2) %>%
  left_join(x = ., y = mod4.nlme.varStruct.A, by = "grpA") %>%
  left_join(x = ., y = mod4.nlme.varStruct.B, by = "grpB") %>%
  mutate(sigma         = mod4.nlme$sigma) %>%
  mutate(StandardError = sigma * varStructA * varStructB) %>%
  mutate(Variance      = StandardError ^ 2)


mod5.nlme.VC <- mod5.nlme$modelStruct$varStruct %>%
  coef(unconstrained = FALSE, allCoef = TRUE) %>%
  enframe(name = "grp", value = "varStruct") %>%
  separate(grp, sep = "[.]", into = c("grpA", "grpB")) %>%
  mutate(sigma         = mod5.nlme$sigma) %>%
  mutate(StandardError = sigma * varStruct) %>%
  mutate(Variance      = StandardError ^ 2)


mod1.lme4.VC <- mod1.lme4 %>%
  tidy(effects = "ran_pars", scales = "vcov") %>%
  separate(term, sep = "__", into = c("term", "grp")) %>%
  filter(group == "Residual")


mod1.glmm.VC <- mod1.glmm %>%
  tidy(effects = "ran_pars", scales = "vcov") %>%
  separate(term, sep = "__", into = c("term", "grp")) %>%
  filter(group == "Residual")


mod2.glmm.VC <- mod2.glmm %>%
  tidy(effects = "ran_pars", scales = "vcov") %>%
  separate(term, sep = "__", into = c("term", "grp")) %>%
  filter(group == "unit", term == "var")


mod3.glmm.VC <- mod3.glmm %>%
  tidy(effects = "ran_pars", scales = "vcov") %>%
  separate(term, sep = "__", into = c("term", "grp")) %>%
  filter(group == "unit", term == "var")


mod5.glmm.VC <- mod5.glmm %>%
  tidy(effects = "ran_pars", scales = "vcov") %>%
  separate(term, sep = "__", into = c("term", "grp")) %>%
  filter(group == "unit", term == "var") %>%
  separate(grp, sep = ":", into = c("grpA", "grpB"))


mod1.somm.VC <- summary(mod1.somm)$varcomp 
mod1.somm.VC <- mod1.somm.VC %>% 
  as_tibble(rownames="grp") %>% 
  mutate(grp = str_replace(grp, "\\..*", "")) %>% 
  filter(grp!="block")


mod2.somm.VC <- summary(mod2.somm)$varcomp
mod2.somm.VC <- mod2.somm.VC %>% 
  as_tibble(rownames="grp") %>% 
  mutate(grp = str_replace(grp, "\\..*", "")) %>% 
  filter(grp!="block") 


mod3.somm.VC <- summary(mod3.somm)$varcomp
mod3.somm.VC <- mod3.somm.VC %>% 
  as_tibble(rownames="grp") %>% 
  mutate(grp = str_replace(grp, "\\..*", "")) %>% 
  filter(grp!="block") 


mod5.somm.VC <- summary(mod5.somm)$varcomp
mod5.somm.VC <- mod5.somm.VC %>% 
  as_tibble(rownames="grp") %>% 
  mutate(grp = str_replace(grp, "\\.yield-yield", "")) %>% 
  filter(grp!="block") 


sasvc <- readr::read_delim("SAS/heterogeneous_error_variance_results_VC.txt", delim="\t") %>% 
  filter(CovParm != "block")

mod1.sas.VC <- sasvc %>%
  filter(mod==1) %>% 
  mutate(Group="homoscedastic") %>% 
  dplyr::select(CovParm, Group, Estimate)

mod2.sas.VC <- sasvc %>%
  filter(mod==2) %>% 
  # mutate(Group = str_remove(Group, "date ")) %>% 
  dplyr::select(CovParm, Group, Estimate)

mod3.sas.VC <- sasvc %>%
  filter(mod==3) %>% 
  # mutate(Group = str_remove(Group, "densf ")) %>% 
  dplyr::select(CovParm, Group, Estimate)

mod4.sas.VC <- sasvc %>%
  filter(mod==4) %>% 
  dplyr::select(CovParm, Estimate) %>% 
  bind_rows(tibble(CovParm=c("EXP date", "EXP densf"), Estimate=c(0.5581,sum(-1.29,-0.507,1.06)))) %>% 
  bind_cols(Group=c("21Aug1990","1","2","4",NA,"28Aug1990","8")) %>% 
  mutate(CovParm = str_remove(CovParm, "EXP "))
mod4.sas.varStruct.A <- mod4.sas.VC %>% filter(CovParm=="date")
mod4.sas.varStruct.B <- mod4.sas.VC %>% filter(CovParm=="densf")
mod4.sas.VC <- expand.grid(mod4.sas.varStruct.A$Group,
                           mod4.sas.varStruct.B$Group,
                           stringsAsFactors = FALSE) %>%
  rename(grpA = Var1, grpB = Var2) %>%
  left_join(x = ., y = mod4.sas.VC %>% dplyr::select(-CovParm) %>% rename(EXPdate=Estimate), 
            by = c("grpA"="Group")) %>%
  left_join(x = ., y = mod4.sas.VC %>% dplyr::select(-CovParm) %>% rename(EXPdensf=Estimate), 
            by = c("grpB"="Group")) %>%
  mutate(sigma    = mod4.sas.VC %>% filter(CovParm=="Residual") %>% pull(Estimate)) %>%
  mutate(Variance = sigma * exp(EXPdate) * exp(EXPdensf)) %>% 
  arrange(grpB, grpA)
  
mod5.sas.VC <- expand.grid(mod4.sas.varStruct.A$Group,
                           mod4.sas.varStruct.B$Group,
                           stringsAsFactors = FALSE) %>% 
  rename(grpA = Var1, grpB = Var2) %>% 
  arrange(grpA, grpB) %>% 
  bind_cols(sasvc %>% filter(mod==5) %>% dplyr::select(CovParm, Estimate)) %>% 
  dplyr::select(CovParm, grpA, grpB, Estimate)











AIC.nlme <- aictab(list(mod1.nlme, mod2.nlme, mod3.nlme, mod4.nlme, mod5.nlme)) %>%
  mutate(Deviance = -2 * Res.LL) # compute deviance manually



AIC.lme4 <- aictab(list(mod1.lme4)) %>% # Mods 2-5 are missing
  mutate(Deviance = -2 * Res.LL) # compute deviance manually


AIC.glmm <- aictab(list(mod1.glmm, mod2.glmm, mod3.glmm, mod5.glmm), 
                   modnames=c("Mod1","Mod2","Mod3","Mod5")) %>% # Mod4 is missing
  mutate(Deviance = -2*LL) # compute deviance manually



somm.mods <- list(mod1.somm, mod2.somm, mod3.somm, mod5.somm)

AIC.somm <- tibble(
    Modnames = paste0("Mod", c(1:3, 5)), # Mod4 is missing
    AIC = somm.mods %>% map("AIC") %>% unlist,
    LL  = somm.mods %>% map("monitor") %>% lapply(. %>% `[`(1, ncol(.))) %>% unlist) %>% # last element in first row of "monitor"
  mutate(Deviance = -2 * LL) %>% # compute deviance manually
  arrange(AIC)


AIC.sas <- readr::read_delim("SAS/heterogeneous_error_variance_results_AIC.txt", delim="\t") %>% 
  mutate(Descr=str_remove(Descr, " \\(kleiner ist besser\\)"),
         Modnames=paste0("Mod",mod)) %>% 
  dplyr::select(-mod) %>% 
  pivot_wider(names_from="Descr", values_from="Value") %>% 
  rename(Deviance = `-2 Res Log-Likelihood`)


plyr::join_all(list(mod1.nlme.VC %>% dplyr::select(grp, Variance) %>% rename(nlme=Variance),
                    mod1.lme4.VC %>% mutate(grp="homoscedastic") %>% 
                      dplyr::select(grp, estimate) %>% rename(lme4=estimate),
                    mod1.glmm.VC %>% mutate(grp="homoscedastic") %>% 
                      dplyr::select(grp, estimate) %>% rename(glmmTMB=estimate),
                    mod1.somm.VC %>% mutate(grp="homoscedastic") %>% 
                      dplyr::select(grp, VarComp) %>% rename(sommer=VarComp),
                    mod1.sas.VC %>% mutate(grp="homoscedastic") %>% 
                      dplyr::select(grp, Estimate) %>% rename(SAS=Estimate)
                    ), 
               by="grp", type="left") %>%
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(mod2.nlme.VC %>% dplyr::select(grp, Variance) %>% rename(nlme=Variance),
                    mod2.nlme.VC %>% dplyr::select(grp) %>% mutate(lme4=NA),
                    mod2.glmm.VC %>% dplyr::select(grp, estimate) %>% 
                      mutate(grp = str_remove(grp, "date")) %>% rename(glmmTMB=estimate),
                    mod2.somm.VC %>% dplyr::select(grp, VarComp) %>% 
                      mutate(grp = str_replace(grp, ":units", "")) %>% rename(sommer=VarComp),
                    mod2.sas.VC %>% 
                      mutate(grp = str_remove(Group, "date ") %>% str_replace("UG","ug")) %>% 
                      dplyr::select(grp, Estimate) %>% rename(SAS=Estimate)
                    ), 
               by="grp", type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(mod3.nlme.VC %>% dplyr::select(grp, Variance) %>% rename(nlme=Variance),
                    mod3.nlme.VC %>% dplyr::select(grp) %>% mutate(lme4=NA),
                    mod3.glmm.VC %>% dplyr::select(grp, estimate) %>% 
                      mutate(grp = str_remove(grp, "densf")) %>% rename(glmmTMB=estimate),
                    mod3.somm.VC %>% dplyr::select(grp, VarComp) %>% 
                      mutate(grp = str_replace(grp, ":units", "")) %>% rename(sommer=VarComp),
                    mod3.sas.VC %>% 
                      mutate(grp = str_remove(Group, "densf ")) %>% 
                      dplyr::select(grp, Estimate) %>% rename(SAS=Estimate)
                    ), 
               by="grp", type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(mod4.nlme.VC %>% dplyr::select(grpA, grpB, Variance) %>% rename(nlme=Variance),
                    mod4.nlme.VC %>% dplyr::select(grpA, grpB) %>% mutate(lme4=NA),
                    mod4.nlme.VC %>% dplyr::select(grpA, grpB) %>% mutate(glmmTMB=NA),
                    mod4.nlme.VC %>% dplyr::select(grpA, grpB) %>% mutate(sommer="in progress"),
                    mod4.sas.VC %>% dplyr::select(grpA, grpB, Variance) %>% rename(SAS=Variance)), 
               by=c("grpA", "grpB"), type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(mod5.nlme.VC %>% dplyr::select(grpA, grpB, Variance) %>% rename(nlme=Variance),
                    mod5.nlme.VC %>% dplyr::select(grpA, grpB) %>% mutate(lme4=NA),
                    mod5.glmm.VC %>% mutate(grpA = str_remove(grpA, "date"),
                                            grpB = str_remove(grpB, "densf")) %>% 
                      dplyr::select(grpA, grpB, estimate) %>% rename(glmmTMB=estimate),
                    mod5.somm.VC %>% mutate(grp = str_replace(grp, ":units", "")) %>% 
                     dplyr::select(grp, VarComp) %>% separate(grp, sep="\\:", into=c("grpA","grpB")) %>%
                     rename(sommer=VarComp),
                    mod5.sas.VC %>% dplyr::select(grpA, grpB, Estimate) %>% rename(SAS=Estimate)), 
               by=c("grpA", "grpB"), type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(AIC.nlme %>% dplyr::select(Modnames, AICc) %>% rename(nlme=AICc),
                    AIC.lme4 %>% dplyr::select(Modnames, AICc) %>% rename(lme4=AICc),
                    AIC.glmm %>% dplyr::select(Modnames, AICc) %>% rename(glmmTMB=AICc),
                    AIC.somm %>% dplyr::select(Modnames, AIC) %>% rename(sommer=AIC),
                    AIC.sas  %>% dplyr::select(Modnames, AICC) %>% rename(sas=AICC)
                    ), by="Modnames", type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)
