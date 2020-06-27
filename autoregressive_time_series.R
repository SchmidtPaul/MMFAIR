params <-
list(hilang = "sas")


# take a character vector of parameters and inject
#   the appropriate script tag for code mirror
# ensures that the script tags are only inserted once
for(i in seq_along(params$hilang)) {
  js_mode <- paste0("\n<script src=\"https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.31.0/mode/", tolower(params$hilang[i]), "/", tolower(params$hilang[i]), ".min.js\"></script>\n")
  cat(htmltools::htmlPreserve(js_mode))
}

knitr::knit_hooks$set(source = function(x, options) {
  if (!is.null(options$hilang)) {
    textarea_id <- paste(sample(LETTERS, 5), collapse = "")
    code_open <- paste0("\n\n<textarea id=\"", textarea_id, "\">\n")
    code_close <- "\n</textarea>"
    jscript_editor <- paste0("\n<script> var codeElement = document.getElementById(\"", textarea_id, "\"); var editor = null; if (null != codeElement) { editor = CodeMirror.fromTextArea(codeElement, { lineNumbers: true, readOnly: true, viewportMargin: Infinity, mode: 'text/x-", tolower(options$hilang), "' }); } </script>\n")
    
    # if the option from_file is set to true then assume that
    #   whatever is in the code chunk is a file path
    if (!is.null(options$from_file) && options$from_file) {
      code_body <- readLines(file.path(x))   
    } else {
      code_body <- x
    }
    
    knitr::asis_output(
      htmltools::htmlPreserve(
        stringr::str_c(
          code_open,
          paste(code_body, collapse = "\n"),
          code_close,
          jscript_editor
        )
      )
    )
  } else {
    stringr::str_c("\n\n```", tolower(options$engine), "\n", paste0(x, collapse = "\n"), "\n```\n\n")
  }
})

# packages
pacman::p_load(dplyr, purrr, tibble, tidyr, stringr, # data handling
               ggplot2, viridis,            # plot
               nlme, lme4, glmmTMB, sommer, # mixed modelling
               AICcmodavg, broom.mixed)     # mixed model extractions

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
  tidy(effects = "ran_pars", scales = "vcov") %>%
  separate(term, sep = "__", into = c("term", "grp"))




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
mod.iid.glmm.VC <- mod.ar1.glmm %>%
  tidy(effects = "ran_pars", scales ="sdcor") %>%
  separate(term, sep = "__", into = c("term", "grp"))
