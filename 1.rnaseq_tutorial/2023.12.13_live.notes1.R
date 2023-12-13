library(tidyverse)
library(limma)
library(kimma)

# data
load("data/dat_tcell.RData")
names(dat_tcell)

# Limma
dat_tcell$targets %>% View()

mm_limma <- model.matrix(~mtb, 
                  data = dat_tcell$targets)
mm_limma

dim(dat_tcell)
fit_limma <- lmFit(object = dat_tcell$E,
                   design = mm_limma,
                   weights = dat_tcell$weights)
efit_limma <- eBayes(fit = fit_limma)
fdr_limma <- topTable(fit = efit_limma)
fdr_limma

fdr_limma <- topTable(fit = efit_limma, 
                      coef = "mtbMtb",
                      number = Inf)
dim(fdr_limma)

fdr_limma <- extract_lmFit(design = mm_limma,
                           fit = efit_limma)
length(fdr_limma)
typeof(fdr_limma)
fdr_limma$lm.fit %>%  head

# Paired design in limma
consensus.corr <- duplicateCorrelation(
  object = dat_tcell$E,
  design = mm_limma,
  block = dat_tcell$targets$ptID
)
consensus.corr$consensus.correlation

# kimma

fit_kimma <- kmFit(dat = dat_tcell,
                   model = "~mtb+(1|ptID)",
                   use_weights = TRUE,
                   run_lm = TRUE, 
                   run_lme = TRUE,
                   metrics = TRUE)
?kmFit
fit_kimma %>% head

fit_kimma_all <- full_join(fit_kimma$lm.fit, 
                           fit_kimma$lme.fit, 
                           by = c("gene"), 
                           suffix = c("_lm","_lme")) %>% 
  #create color variable
  mutate(diff = AIC_lme-AIC_lm,
         diff_col = case_when(diff<=-7 | diff>=7 ~ "Strong",
                              diff<=-2 | diff>=2 ~ "Moderate",
                              TRUE~"No difference"),
         diff_col = factor(diff_col, 
                           levels=c("No difference",
                                    "Moderate","Strong")))
#Comparing model fits
fit_kimma_all %>%
  ggplot(aes(x = AIC_lm, y = AIC_lme)) +
  geom_point(alpha = 0.2, aes(color=diff_col)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  coord_fixed() +
  labs(title = "AIC", color="AIC difference") +
  scale_color_manual(values=c("grey40","orange","darkred")) +
  annotate("text", x=150, y=0, label="Better fit by lme")+
  annotate("text", x=0, y=150, label="Better fit by lm")

mean(fit_kimma$lm.fit$AIC)
mean(fit_kimma$lme.fit$AIC)

# Significant genes per model
summarise_kmFit(fdr_limma$lm)
summarise_kmFit(fit_kimma$lm)
summarise_kmFit(fit_kimma$lme)
