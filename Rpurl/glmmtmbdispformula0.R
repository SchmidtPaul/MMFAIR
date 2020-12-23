pacman::p_load(conflicted,
               tidyverse, 
               nlme, glmmTMB,
               broom.mixed,
               emo, flair)

# package function conflicts
conflict_prefer("filter", "dplyr")











dat <- agridat::mcconway.turnip %>% 
  mutate(unit = 1:n()) %>% 
  mutate_at(vars(density, unit), as.factor)

StandErrMod <- glmmTMB(
  yield ~ 
    gen*date*density +
    (1 | block),
  

  REML = TRUE,
  data = dat
)

PseudErrMod <- glmmTMB(
  yield ~
    gen*date*density +
    (1 | block) +
    (1 | unit),      # Pseudo Err
  dispformula = ~ 0, # ErrVar = 0
  REML = TRUE,
  data = dat
)

## mixedup::extract_vc(StandErrMod)


## mixedup::extract_vc(PseudErrMod)


## AICcmodavg::aictab(list(StandErrMod, PseudErrMod),
##                    c("StandErrMod", "PseudErrMod"),
##                    second.ord = FALSE)


diag_glmmTMB <- glmmTMB(
  yield ~ 
    gen*date*density +
    (1 | block) +
    diag(date + 0 | unit),
  dispformula = ~ 0,
  REML = TRUE,
  data = dat
)

diag_lme <- lme(
  yield ~ 
    gen*date*density, 
  random  = ~ 1 | block,
  weights = varIdent(form =  ~ 1 | date),
  
  
  data    = dat
)

glmmTMB_vc <- bind_rows(

diag_glmmTMB %>%
  tidy(effects = "ran_pars", scales = "vcov") %>%
  filter(group == "block") %>% 
  mutate(grp = str_remove(term, "var__")) %>% 
  select(group, grp, estimate) %>% 
  rename(variance = estimate,
         effect = group)

,

diag_glmmTMB %>% 
  mixedup::extract_cor_structure(which_cor="diag") %>% 
  pivot_longer( cols = 2:3, values_to ="variance", names_to="grp") %>% 
  mutate(variance = variance ^ 2) %>% 
  rename(effect = group)

,

tibble(effect = "Residual",
       grp = NA_character_,
       variance = glance(diag_glmmTMB) %>% pull(sigma) %>% `^`(2))

)



lme_vc <- bind_rows(

diag_lme %>% 
  tidy(effects = "ran_pars", scales = "vcov") %>% 
  filter(group=="block") %>% 
  mutate(grp = str_remove(term, "var_")) %>% 
  select(group, grp, estimate) %>% 
  rename(variance = estimate,
         effect = group)

,

diag_lme$modelStruct$varStruct %>%
  coef(unconstrained = FALSE, allCoef = TRUE) %>%
  enframe(name = "grp", value = "varStruct") %>%
  mutate(sigma         = diag_lme$sigma) %>%
  mutate(StandardError = sigma * varStruct) %>%
  mutate(variance      = StandardError ^ 2) %>% 
  mutate(effect = "Residual") %>% 
  select(effect, grp, variance)

) 



glmmTMB_fit <- diag_glmmTMB %>% 
  glance() %>% select(logLik:BIC)



lme_fit <- diag_lme %>%   
  glance() %>% select(logLik:BIC) 
