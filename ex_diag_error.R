# packages for better formatting tables for html output
library(kableExtra)
library(formattable)
library(stringr); library(stringi)

# packages
pacman::p_load(dplyr, purrr, tibble, tidyr, # data handling
               nlme, lme4, glmmTMB, sommer, # mixed modelling
               AICcmodavg, broom.mixed)     # mixed model extractions
# data
dat <- agridat::mcconway.turnip %>%
  as_tibble() %>% 
  mutate(densf = density %>% as.factor)

boxplot(yield ~ date,    data=dat)
boxplot(yield ~ density, data=dat)

tibble(Model = paste0("mod", 1:5) %>% cell_spec(bold=T),
       Block    = "Identity",
       Genotype = "Identity",
       Date     = c("Identity", "Identity", "Diagonal", "Diagonal", "Diag-"),
       Density  = c("Identity", "Diagonal", "Identity", "Diagonal", "onal"),
       `variance parameters` = c(1,2,4,5,8),
       `variance estimates`  = c(1,2,4,8,8)) %>% 
  mutate(Date    = ifelse(Date    %in% c("Diagonal", "Diag-"), cell_spec(Date,    bold=T), Date),
         Density = ifelse(Density %in% c("Diagonal", "onal"),  cell_spec(Density, bold=T), Density)) %>%
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
                full_width = FALSE) %>% 
  add_header_above(c(" ", "Term in multiplicative variance structure"=4, "Number of error"=2))

dat <- dat %>% 
  mutate(date_densf = interaction(date, densf)) # needed for mod5

mod1.nlme <- nlme::lme(fixed   = yield ~ gen * date * densf, 
                       random  = ~ 1|block,
                       weights = NULL, # default, i.e. homoscedastic errors
                       data    = dat)

mod2.nlme <- mod1.nlme %>% 
  update(weights = varIdent(form=~1|date))

mod3.nlme <- mod1.nlme %>% 
  update(weights = varIdent(form=~1|densf))

mod4.nlme <- mod1.nlme %>% 
  update(weights = varComb(varIdent(form=~1|date), 
                           varIdent(form=~1|densf))) 

mod5.nlme <- mod1.nlme %>% 
  update(weights = varIdent(form=~1|date_densf))

dat <- dat %>% 
  mutate(unit = 1:n() %>% as.factor) # new column with running number

mod1.glmm <- glmmTMB(formula = yield ~ gen * date * densf + (1|block),
                     REML    = TRUE,
                     data    = dat) # default = homoscedastic error variance

mod1b.glmm <- mod1.glmm %>% 
  update(.~. + (1|unit), # add random term to mimic homoscedastic error variance
         dispformula=~0) # fix original error variance to 0

mod2.glmm <- mod1.glmm %>% 
  update(.~. + diag(date+0|unit), dispformula=~0)

mod3.glmm <- mod1.glmm %>% 
  update(.~. + diag(densf+0|unit), dispformula=~0)

# mod4.glmm ?

mod5.glmm <- mod1.glmm %>% 
  update(.~. + diag(date:densf+0|unit), dispformula=~0)

mod1.glmm %>% VarCorr()

mod1b.glmm %>% VarCorr()

mod1.nlme.VC <- tibble(grp="homoscedastic", varStruct=1) %>% 
  mutate(sigma         = mod1.nlme$sigma) %>% 
  mutate(StandardError = sigma*varStruct) %>% 
  mutate(Variance      = StandardError^2)
mod1.nlme.VC %>%
  mutate_if(is.double, round, 3) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

mod2.nlme.VC <- mod2.nlme$modelStruct$varStruct %>% 
  coef(unconstrained=FALSE, allCoef=TRUE) %>% 
  enframe(name="grp", value="varStruct") %>% 
  mutate(sigma         = mod2.nlme$sigma) %>% 
  mutate(StandardError = sigma*varStruct) %>% 
  mutate(Variance      = StandardError^2)
mod2.nlme.VC %>%
  mutate_if(is.double, round, 3) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

mod3.nlme.VC <- mod3.nlme$modelStruct$varStruct %>% 
  coef(unconstrained=FALSE, allCoef=TRUE) %>% 
  enframe(name="grp", value="varStruct") %>% 
  mutate(sigma         = mod3.nlme$sigma) %>% 
  mutate(StandardError = sigma*varStruct) %>% 
  mutate(Variance      = StandardError^2)
mod3.nlme.VC %>%
  mutate_if(is.double, round, 3) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

mod4.nlme.varStruct.A <- mod4.nlme$modelStruct$varStruct$A %>% 
  coef(unconstrained=FALSE, allCoef=TRUE) %>% 
  enframe(name="grpA", value="varStructA")

mod4.nlme.varStruct.B <- mod4.nlme$modelStruct$varStruct$B %>% 
  coef(unconstrained=FALSE, allCoef=TRUE) %>% 
  enframe(name="grpB", value="varStructB")

mod4.nlme.VC <- expand.grid(mod4.nlme.varStruct.A$grpA, 
                            mod4.nlme.varStruct.B$grpB,
                            stringsAsFactors=FALSE) %>% 
  rename(grpA=Var1, grpB=Var2) %>% 
  left_join(x=., y=mod4.nlme.varStruct.A, by="grpA") %>% 
  left_join(x=., y=mod4.nlme.varStruct.B, by="grpB") %>% 
  mutate(sigma         = mod4.nlme$sigma) %>% 
  mutate(StandardError = sigma*varStructA*varStructB) %>% 
  mutate(Variance      = StandardError^2)
mod4.nlme.VC %>%
  arrange(grpA, grpB) %>% 
  mutate_if(is.double, round, 3) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

mod5.nlme.VC <- mod5.nlme$modelStruct$varStruct %>% 
  coef(unconstrained=FALSE, allCoef=TRUE) %>% 
  enframe(name="grp", value="varStruct") %>% 
  separate(grp, sep="[.]", into=c("grpA","grpB")) %>% 
  mutate(sigma         = mod5.nlme$sigma) %>% 
  mutate(StandardError = sigma*varStruct) %>% 
  mutate(Variance      = StandardError^2)
mod5.nlme.VC %>%
  mutate_if(is.double, round, 3) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

mod1.glmm.VC <- mod1.glmm %>% 
  tidy(effects="ran_pars", scales="vcov") %>% 
  separate(term, sep="__", into=c("term","grp")) %>% 
  filter(group=="Residual")
mod1.glmm.VC %>%
  mutate_if(is.double, round, 3) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

mod2.glmm.VC <- mod2.glmm %>% 
  tidy(effects="ran_pars", scales="vcov") %>% 
  separate(term, sep="__", into=c("term","grp")) %>% 
  filter(group=="unit", term=="var")
mod2.glmm.VC %>%
  mutate_if(is.double, round, 3) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

mod3.glmm.VC <- mod3.glmm %>% 
  tidy(effects="ran_pars", scales="vcov") %>% 
  separate(term, sep="__", into=c("term","grp")) %>% 
  filter(group=="unit", term=="var")
mod3.glmm.VC %>%
  mutate_if(is.double, round, 3) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

mod5.glmm.VC <- mod5.glmm %>% 
  tidy(effects="ran_pars", scales="vcov") %>% 
  separate(term, sep="__", into=c("term","grp")) %>% 
  filter(group=="unit", term=="var") %>% 
  separate(grp, sep=":", into=c("grpA","grpB"))
mod5.glmm.VC %>%
  mutate_if(is.double, round, 3) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

mod1.somm <- mmer(fixed  = yield ~ gen*date*densf,
                  random = ~ block,
                  rcov   = ~ units, # default
                  data   = dat)

mod2.somm <- mmer(fixed  = yield ~ gen*date*densf,
                  random = ~ block,
                  rcov   = ~ vs(ds(date),units),
                  data   = dat)

mod3.somm <- mmer(fixed  = yield ~ gen*date*densf,
                  random = ~ block,
                  rcov   = ~ vs(ds(densf),units),
                  data   = dat)

mod4.somm <- mmer(fixed  = yield ~ gen*date*densf,
                  random = ~ block,
                  rcov   = ~ vs(ds(date),ds(densf),units),
                  data   = dat)

mod5.somm <- mmer(fixed  = yield ~ gen*date*densf,
                  random = ~ block,
                  rcov   = ~ vs(ds(date_densf),units),
                  data   = dat)

mod1.somm.VC <- summary(mod1.somm)$varcomp %>% 
  as_tibble(rownames="grp") %>% 
  mutate(grp = str_replace(grp, "\\..*", "")) %>% 
  filter(grp!="block")

mod2.somm.VC <- summary(mod2.somm)$varcomp %>% 
  as_tibble(rownames="grp") %>% 
  mutate(grp = str_replace(grp, "\\..*", "")) %>% 
  filter(grp!="block")

mod3.somm.VC <- summary(mod3.somm)$varcomp %>% 
  as_tibble(rownames="grp") %>% 
  mutate(grp = str_replace(grp, "\\..*", "")) %>% 
  filter(grp!="block")

mod4.somm.VC <- summary(mod4.somm)$varcomp %>% 
  as_tibble(rownames="grp") %>% 
  mutate(grp = str_replace(grp, "\\..*", "")) %>% 
  filter(grp!="block")

mod5.somm.VC <- summary(mod5.somm)$varcomp %>% 
  as_tibble(rownames="grp") %>% 
  mutate(grp = str_replace(grp, "\\.yield-yield", "")) %>% 
  filter(grp!="block")

row.names(summary(mod2.somm)$varcomp)

AIC.nlme <- aictab(list(mod1.nlme, mod2.nlme, mod3.nlme, mod4.nlme, mod5.nlme)) %>%
  mutate(Deviance = -2*Res.LL) # compute deviance

AIC.nlme %>%  
  dplyr::select(Modnames, K, AICc, Delta_AICc, Res.LL, Deviance) %>% 
  mutate_at(vars(AICc:Deviance), round, 1) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

AIC.glmm <- aictab(list(mod1.glmm, mod2.glmm, mod3.glmm, mod5.glmm), 
                   modnames=c("Mod1","Mod2","Mod3","Mod5")) %>% # Mod4 is missing
  mutate(Deviance = -2*LL) # compute deviance

AIC.glmm %>%  
  dplyr::select(Modnames, K, AICc, Delta_AICc, LL, Deviance) %>%  
  mutate_at(vars(AICc:Deviance), ~round(., 1)) %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(mod1.nlme.VC %>% dplyr::select(grp, Variance) %>% rename(nlme=Variance),
                    mod1.nlme.VC %>% dplyr::select(grp) %>% mutate(lme4=NA),
                    mod1.glmm.VC %>% mutate(grp="homoscedastic") %>% 
                      dplyr::select(grp, estimate) %>% rename(glmmTMB=estimate),
                    mod1.nlme.VC %>% dplyr::select(grp) %>% mutate(sommer="in progress")), 
               by="grp", type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(mod2.nlme.VC %>% dplyr::select(grp, Variance) %>% rename(nlme=Variance),
                    mod2.nlme.VC %>% dplyr::select(grp) %>% mutate(lme4=NA),
                    mod2.glmm.VC %>% dplyr::select(grp, estimate) %>% 
                      mutate(grp = str_remove(grp, "date")) %>% rename(glmmTMB=estimate),
                    mod2.nlme.VC %>% dplyr::select(grp) %>% mutate(sommer="in progress")), 
               by="grp", type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(mod3.nlme.VC %>% dplyr::select(grp, Variance) %>% rename(nlme=Variance),
                    mod3.nlme.VC %>% dplyr::select(grp) %>% mutate(lme4=NA),
                    mod3.glmm.VC %>% dplyr::select(grp, estimate) %>% 
                      mutate(grp = str_remove(grp, "densf")) %>% rename(glmmTMB=estimate),
                    mod3.nlme.VC %>% dplyr::select(grp) %>% mutate(sommer="in progress")), 
               by="grp", type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(mod4.nlme.VC %>% dplyr::select(grpA, grpB, Variance) %>% rename(nlme=Variance),
                    mod4.nlme.VC %>% dplyr::select(grpA, grpB) %>% mutate(lme4=NA),
                    mod4.nlme.VC %>% dplyr::select(grpA, grpB) %>% mutate(glmmTMB=NA),
                    mod4.nlme.VC %>% dplyr::select(grpA, grpB) %>% mutate(sommer="in progress")), 
               by=c("grpA", "grpB"), type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(mod5.nlme.VC %>% dplyr::select(grpA, grpB, Variance) %>% rename(nlme=Variance),
                    mod5.nlme.VC %>% dplyr::select(grpA, grpB) %>% mutate(lme4=NA),
                    mod5.glmm.VC %>% mutate(grpA = str_remove(grpA, "date"),
                                            grpB = str_remove(grpB, "densf")) %>% 
                      dplyr::select(grpA, grpB, estimate) %>% rename(glmm=estimate),
                   mod5.nlme.VC %>% dplyr::select(grpA, grpB) %>% mutate(sommer="in progress")), 
               by=c("grpA", "grpB"), type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)

plyr::join_all(list(AIC.nlme %>% dplyr::select(Modnames, AICc) %>% rename(nlme=AICc),
                    AIC.nlme %>% dplyr::select(Modnames) %>% mutate(lme4=NA),
                    AIC.glmm %>% dplyr::select(Modnames, AICc) %>% rename(glmm=AICc),
                    AIC.nlme %>% dplyr::select(Modnames) %>% mutate(lme4="in progress")
                    ), by="Modnames", type="left") %>% 
  kable(escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("bordered", "hover", "condensed", "responsive"), 
                full_width = FALSE)
