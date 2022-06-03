## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2022-03-22
##
## Copyright (c) Jan D. Lanzer, 2022
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## Evaluate Plasma levels with clinical parameters
## nur mit hf patients

library(ComplexHeatmap)
library(tidyverse)
library(haven)
library(tidyverse)

meta = read_sav(file = "data/plasma/2200301_Data_AFFECT.sav")
anp= read.csv("data/plasma/ANGPTL4_20220316.csv")
meta%>%select(starts_with("NYHA"))
hist(meta$NYHA_change)

#map= read.csv("data/plasma/pair_matching.csv")

df= meta %>%
  inner_join(anp%>% rename(Pat_ID= Sample_clean)) %>%
  mutate(hfpef = factor(hfpef, levels=c(0,1)))%>%
  as_factor()%>%
  rename(anptl4= unverd..Probe)



# ## 1. use numerical values for correlation testing ------------------------------------------

## Test correlation for all patients

colnames(df)
nums= df%>% dplyr::select(where(is.numeric))
#remove kidney fail stage because of number of NAs
cols.r= c("kidney_fail_stage")
test.cols= colnames(nums)
test.cols=test.cols[!test.cols %in% cols.r]

cor.df= sapply(test.cols, function(x){
  print(x)
  res= cor.test(df[[x]], df$anptl4, method = "spearman")
  c("p"= res$p.value, res$estimate)
})%>% t()

cor.df= cor.df[!rowSums(!is.finite(cor.df)),]

cor.sig= cor.df %>%as.data.frame() %>% rownames_to_column("feat")%>% as_tibble() %>% arrange(p)
cor.sig%>%
  filter(p<0.1)%>% arrange(p)%>% print(n=100)

## Test correlation for HF patients
df_hf= df %>% filter(hfpef ==1)

cor.df_hf= sapply(test.cols, function(x){
  print(x)
  #s= shapiro.test(df_hf[[x]])$p.value
  #print(paste0("pval_normal", s>0.1))
  res= cor.test(df_hf[[x]], df_hf$anptl4,method = "pearson")
  c("p"= res$p.value, res$estimate)
})%>% t()

cor.df_hf= cor.df_hf[!rowSums(!is.finite(cor.df_hf)),]

cor.sig_hf= cor.df_hf %>%as.data.frame() %>% rownames_to_column("feat")%>% as_tibble()
cor.sig_hf %>% filter(grepl("E_e", feat))
cor.sig_hf%>%
  filter(p<0.1)%>% arrange(p)%>% print(n=100)

ggplot(df_hf, aes(x= anptl4, y= Strain_Baseline))+
  geom_point()
cor.test(df_hf$anptl4, df_hf$Strain_Baseline, method = "spearman")


## 2. Wilcox for categorical 2 level
fact= df%>% dplyr::select(where(is.factor))
fact.l= sapply(fact, function(x) (length(unique(x))))
test.cols= names(fact.l[fact.l <3])

t.df= sapply(test.cols, function(x){
  print(x)
  v= unique(df[[x]])
  x1= df%>% filter(!!as.symbol(x) == v[1])%>% pull(anptl4)
  x2= df%>% filter(!!as.symbol(x)== v[2])%>% pull(anptl4)
  if(length(x1>3) & length(x2)>4){
    res= wilcox.test(x1, x2, paired = F)
    return(  c("p"= res$p.value, res$statistic))
  }else{ return(rep(NA, 2))}
})%>% t()

t.df= t.df[!rowSums(!is.finite(t.df)),]

t.df %>%as.data.frame() %>% rownames_to_column("feat")%>% as_tibble() %>% arrange(p)%>%
  print(n=199)

cor.sig= t.df %>%as.data.frame() %>% rownames_to_column("feat")%>% as_tibble() %>%
  #filter(p<0.1)%>%
  arrange(p)%>% print(n=100)


## 3. ANOVA for categeorical >2 level
test.cols= names(fact.l[fact.l >2])
test.cols = c(test.cols, "NYHA_stage_base", "NYHA_change", "NYHA_stage_12MFU")
t.df= sapply(test.cols, function(x){
  print(x)
  res= aov(formula=as.formula(paste0("anptl4 ~ ", x)), data=df)
  sum_test = unlist(summary(res))


})%>% t()
t.df


#save nyha data for plotting.


df %>% select(starts_with("NYHA"), Pat_ID, hfpef ,anptl4,logntprobnp )%>% write.csv(., "output/nyha_patient_data.csv")
plot(df$NYHA_stage_base, df$NYHA_change)
## 4.LR
library(pROC)

lr.bnp= glm(hfpef~ logntprobnp, data = df, family= "binomial")
auc(df$hfpef, lr.bnp$fitted.values)

lr.bnp.ang= glm(hfpef~ anptl4, data = df, family= "binomial")
auc(df$hfpef, lr.bnp.ang$fitted.values)

lm()
lr.bnp.ang= glm(hfpef~ logntprobnp+age+sex+anptl4, data = df, family= "binomial")
coef(summary(lr.bnp.ang))
coef(summary(lr.bnp))
coefficients(lr.bnp)
lr.bnp$fitted.values
lr.bnp$aic
coef(summary(lr.bnp))


df2 = df2 %>% mutate(anptl4 = scale(anptl4),
                     logntprobnp = scale(logntprobnp))

df2= df2 %>% mutate(nyha.m= factor(ifelse(NYHA_stage_base>1, 1, 0)))

lr.bnp.ang= glm(nyha.m~logntprobnp+ anptl4, data = df2, family= "binomial")
auc(df2$nyha.m, lr.bnp.ang$fitted.values)
coef(summary(lr.bnp.ang))

df2= df %>% mutate(nyha.m= factor(ifelse(NYHA_change>0, 1, 0)))

lr.bnp.ang= glm(nyha.m~ logntprobnp, data = df2, family= "binomial")
auc(df2$nyha.m, lr.bnp.ang$fitted.values)
coef(summary(lr.bnp.ang))

lr.bnp.ang= glm(nyha.m~ logntprobnp+anptl4, data = df2, family= "binomial")
coef(summary(lr.bnp.ang))
lr.bnp= glm(nyha.m~ logntprobnp, data = df2, family= "binomial")
coef(summary(lr.bnp))
lr.ang= glm(nyha.m~ anptl4, data = df2, family= "binomial")
coef(summary(lr.ang))
auc(df2$nyha.m, lr.ang$fitted.values)
auc(df2$nyha.m, lr.bnp$fitted.values)

coef(summary(lr.bnp.ang))

library(MASS)
model <- lda(NYHA_stage_base~logntprobnp+anptl4, data = df2)
model
coef(summary(model))
plot(model)
x= apply(model$posterior, 1, max)
auc(df2$NYHA_stage_base ,x)


# Add plots -----------------------------------------------------------------------------------
## significant correlations in hf patients
cor.feat= cor.sig_hf%>% filter(p<0.05)%>% pull(feat)
pls= map(cor.feat[1:6], function(x){

  lm.fit= lm(paste0( "anptl4 ~", x), data= df_hf)
  r2= round(summary(lm.fit)$r.squared, 2)
  p.= round(summary(lm.fit)$coefficient[2, 4],2)

  ggplot(df_hf, aes(y= anptl4, x = !!as.symbol(x))) +
    geom_point() +
    stat_smooth(method = "lm", col = "red")+
    geom_text(label = paste0("R²=", r2 , "\n", "p-val=", p.), x = min(df_hf[[x]]) , y= 150, size= 3)


})
names(pls)= cor.feat[1:6]
cowplot::plot_grid(plotlist= pls)

pls= map(c("Strain_Baseline", "TAPSE_Baseline", "sves_6mfu"), function(x){
  print(x)
  lm.fit= lm(paste0( "anptl4 ~", x), data= df_hf)
  r2= round(summary(lm.fit)$r.squared, 2)
  p.= round(summary(lm.fit)$coefficient[2, 4],2)
  print(c(r2, p.))
  print(min(df_hf[[x]]))
  ggplot(df_hf, aes(y= anptl4, x = !!as.symbol(x))) +
    geom_point() +
    stat_smooth(method = "lm", col = "red")+
    geom_text(label = paste0("R²=", r2 , "\n", "p-val=", p.),
              x = (min(df_hf[[x]], na.rm = T))+0.1* min(df_hf[[x]], na.rm = T) ,
              y= 155, size= 3)


})

x= "Strain_Baseline"
  print(x)
  lm.fit= lm(paste0( "anptl4 ~", x), data= df_hf)
  r2= round(summary(lm.fit)$r.squared, 2)
  p.= round(summary(lm.fit)$coefficient[2, 4],2)
  print(c(r2, p.))
  print(min(df_hf[[x]]))
  p1 = ggplot(df_hf, aes(y= anptl4, x = !!as.symbol(x))) +
    geom_point() +
    stat_smooth(method = "lm", col = "red")+
    geom_text(label = paste0("R²=", r2 , "\n", "p-val=", p.),
              x = (min(df_hf[[x]], na.rm = T))+0.1* min(df_hf[[x]], na.rm = T) ,
              y= 155, size= 3)
x= "TAPSE_Baseline"
lm.fit= lm(paste0( "anptl4 ~", x), data= df_hf)
r2= round(summary(lm.fit)$r.squared, 2)
p.= round(summary(lm.fit)$coefficient[2, 4],2)
print(c(r2, p.))
print(min(df_hf[[x]]))
p2 = ggplot(df_hf, aes(y= anptl4, x = !!as.symbol(x))) +
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  geom_text(label = paste0("R²=", r2 , "\n", "p-val=", p.),
            x = (min(df_hf[[x]], na.rm = T))+0.2* min(df_hf[[x]], na.rm = T) ,
            y= 155, size= 3)


p2
x= "sves_6mfu"
lm.fit= lm(paste0( "anptl4 ~", x), data= df_hf)
r2= round(summary(lm.fit)$r.squared, 2)
p.= round(summary(lm.fit)$coefficient[2, 4],2)
print(c(r2, p.))
print(min(df_hf[[x]]))
p3 = ggplot(df_hf, aes(y= anptl4, x = !!as.symbol(x))) +
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  geom_text(label = paste0("R²=", r2 , "\n", "p-val=", p.),
            x = 1000,
            y= 50, size= 3)+
  labs(y= "ANGTPL4 pg/ml")

p3
p.c= cowplot::plot_grid(p1+labs(y= "ANGTPL4 pg/ml"),
                        p2+labs(y= "ANGTPL4 pg/ml"),
                        p3+labs(y= "ANGTPL4 pg/ml"), nrow = 1)
pdf("output/figures/main/angptl4_clinvar.pdf",
    height= 3,
    width= 8)
p.c
dev.off()


ggplot(df, aes(x= hfpef, y=LVEF_base))+
  geom_boxplot()+
  geom_jitter()

ggplot(df, aes(x= factor(NYHA_stage_base), y=LVEF_base))+
  geom_boxplot()+
  geom_jitter()

ggplot(df, aes(NYHA_stage_base, ..count..)) +
  geom_bar(aes(fill = hfpef), position = "dodge")

ggplot(df, aes(x= factor(NYHA_stage_base), y=hfpef))+
  geom_boxplot()+
  geom_jitter()

ggplot(df, aes(x= hfpef, y= anptl4))+
  geom_boxplot()+
  geom_jitter()

ggplot(df, aes(x= dyspnoe , y= anptl4))+
  geom_boxplot()+
  geom_jitter()


ggplot(df, aes(x= hfpef, y= logntprobnp))+
  geom_boxplot()
df$LVEF_base
ggplot(df, aes(x= Strain_Baseline, y= LVEF_base))+
  geom_point()
ggplot(df_hf, aes(x= Strain_Baseline, y= LVEF_base))+
  geom_point()
ggplot(df, aes(x= anptl4, y= BMI))+
  geom_point()

ggplot(df, aes(x= diabetes, y= anptl4))+
  geom_boxplot()

ggplot(df, aes(x= E_A_Ratio_Baseline, y= anptl4, color = hfpef))+
  geom_point()
ggplot(df, aes(x= E_A_Ratio_Baseline, y= anptl4, color = hfpef))+
  geom_point()

ggplot(df, aes(x= E_e_ratio_12MFU, y= anptl4))+
  geom_point()

df$E_e_ratio_12MFU
df$E_e_peak_12MFU

df= df %>% mutate(bnp.diff= ntprobnp_12MFU-ntprobnp)

ggplot(df, aes(x= anptl4, y= bnp.diff))+
  geom_point()
cor.test(df$bnp.diff, df$anptl4)

p_NYHA_anptl4= cowplot::plot_grid(
ggplot(df, aes(x= factor(NYHA_stage_base), y= anptl4))+
  geom_boxplot()+
  geom_jitter()+
  labs(x= "NYHA baseline",
       y= "Angptl4 pg/ml")+theme_bw(),
ggplot(df, aes(x= factor(NYHA_stage_12MFU), y= anptl4))+
  geom_boxplot()+
  geom_jitter()+
  labs(x= "NYHA 12months",
       y= "Angptl4 pg/ml")+
  theme_bw()
)
p_NYHA_anptl4_HF= cowplot::plot_grid(
  ggplot(df_hf, aes(x= factor(NYHA_stage_base), y= anptl4 ))+
    geom_boxplot()+
    geom_jitter()+
    labs(x= "NYHA baseline",
         y= "Angptl4 pg/ml")+theme_bw(),
  ggplot(df_hf, aes(x= factor(NYHA_stage_12MFU), y= anptl4))+
    geom_boxplot()+
    geom_jitter()+
    labs(x= "NYHA 12months",
         y= "Angptl4 pg/ml")+
    theme_bw()
)


p_NYHA_anptl4

cowplot::plot_grid(
  ggplot(df_hf, aes(x= factor(NYHA_stage_base), y= logntprobnp))+
    geom_boxplot()+
    geom_jitter(),
  ggplot(df_hf, aes(x= factor(NYHA_change), y= logntprobnp))+
    geom_boxplot()+
    geom_jitter(),
  ggplot(df_hf, aes(x= factor(NYHA_stage_12MFU), y= logntprobnp))+
    geom_boxplot()+
    geom_jitter()
)

cowplot::plot_grid(
  ggplot(df_hf, aes(x= factor(NYHA_stage_base), y= anptl4))+
    geom_boxplot()+
    geom_jitter(),
  ggplot(df_hf, aes(x= factor(NYHA_change), y= anptl4))+
    geom_boxplot()+
    geom_jitter(),
  ggplot(df_hf, aes(x= factor(NYHA_stage_12MFU), y= anptl4))+
    geom_boxplot()+
    geom_jitter()
)

ggplot(df%>%filter(hfpef== 1), aes(x= anptl4, y= Strain_Baseline))+
  geom_point()
df$E_e_ratio_peak_baseline
df$MAPSE_Baseline
ggplot(df, aes(x= anptl4, y= age))+
  geom_point()

ggplot(df, aes(x= factor(EHRA_stage), y= anptl4))+
  geom_boxplot()+
  geom_jitter()
cor.test(df$bnp.diff, df$anptl4)
df$sves_6mfu

