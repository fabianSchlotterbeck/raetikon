library(ggplot2)
library (lme4)
library(lmerTest)
library(Hmisc)
library(blme)

rm(list=ls())

load ("preprocessed_data_panuelerkopf.RData")

plots_dir<- file.path(getwd(),"plots")
if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

#TODO: Maybe move to preprocessing
items_comp_scope <- subset(items_comp_scope, roi == "rt.14")

items_comp_scope$type <- as.factor(ifelse(items_comp_scope$conditions %in% 9:12, "exactly","less"))
contrasts(items_comp_scope$type) <- rbind(-.5,.5)
colnames(contrasts(items_comp_scope$type)) <- levels(items_comp_scope$type)[2]


items_comp_scope$modal <- as.factor(ifelse(items_comp_scope$conditions %in% c(9,10,13,14), "present","absent"))
contrasts(items_comp_scope$modal) <- rbind(-.5,.5)
colnames(contrasts(items_comp_scope$modal)) <- levels(items_comp_scope$modal)[2]

items_comp_scope$question <- as.factor(ifelse(items_comp_scope$conditions %in% c(9,11,13,15), "match","mismatch"))
contrasts(items_comp_scope$question) <- rbind(-.5,.5)
colnames(contrasts(items_comp_scope$question)) <- levels(items_comp_scope$question)[2]

items_comp_scope$response_num <- ifelse(items_comp_scope$response == "True", 1, 0)

ses<-aggregate(response_num~question+modal+type, data=items_comp_scope,
               FUN=function(x) {binconf(x=sum(x), n=length(x))})
agg<-cbind(cbind(ses[1:3],ses[,4][,1],ses[,4][,2],ses[,4][,3]))

items_comp_scope$corr_resp  <- ifelse(items_comp_scope$conditions %in% c(9,11,13,15), 1, 0)

items_comp_scope$acc <- items_comp_scope$corr_resp == items_comp_scope$response_num

aggregate(1-acc~question+modal+type, data=items_comp_scope, FUN=sum)
aggregate(1-acc~modal+type, data=items_comp_scope, FUN=mean)
aggregate(acc~question+modal, data=items_comp_scope, FUN=mean)

names(agg)<-c("question","modal","type", "mean_judgment", "lower", "upper")

eb <- ggplot(data = agg, aes(question, mean_judgment, fill = modal)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9))+
  facet_wrap(~type)+
  ylab("rel. freq. \"yes\" (and approx. 95% CIs)")+
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size=16))
ggsave(file.path(plots_dir,"mean_judgments_comp_scope.pdf"), device="pdf", width=8, height=5)

m.scope.0.0 <- glmer(response_num~question*modal*type
                     +(1|item)+(1|id),
                     family=binomial, data=items_comp_scope,
                     glmerControl(optimizer = c("bobyqa")))
summary(m.scope.0.0)

m.scope.0.1 <- update(m.scope.0.0, .~.-(1|item))
anova(m.scope.0.1,m.scope.0.0)

m.scope.0.2 <- update(m.scope.0.0, .~.-(1|id))
anova(m.scope.0.2,m.scope.0.0)

m.scope.0.3 <- update(m.scope.0.1, .~.+(0+modal|id))
anova(m.scope.0.3,m.scope.0.1)

m.scope.0.4 <- update(m.scope.0.1, .~.+(0+type|id))
anova(m.scope.0.1,m.scope.0.4)


m.scope.0.5 <- update(m.scope.0.1, .~.+(0+question|id))
anova(m.scope.0.5,m.scope.0.1)
summary(m.scope.0.1)

#simplify fixed effects
m.scope.1.1 <- update(m.scope.0.1, .~.-question:modal:type)
anova(m.scope.1.1,m.scope.0.1)
summary(m.scope.0.1)

items_comp_scope_exactly <- subset(items_comp_scope, type=="exactly")

m.scope.exactly.0 <- glmer(response_num~question*modal
                           #+ (0 + modal | id) + (0 + modal | item)
                           #+ (1|item)
                           + (1|id),                         ,
                           family=binomial,
                           data=items_comp_scope_exactly)

summary(m.scope.exactly.0)

m.scope.exactly.1 <- update(m.scope.exactly.0, .~.-question:modal)
anova(m.scope.exactly.1,m.scope.exactly.0)

items_comp_scope_exactly_mismatch <- subset(items_comp_scope_exactly, question=="mismatch")

m.scope.exactly.mismatch.0 <- glmer(response_num~modal+
                                      (1|id),
                                    family=binomial,
                                    glmerControl(optimizer = c("bobyqa")),
                                    data=items_comp_scope_exactly_mismatch)
summary(m.scope.exactly.mismatch.0)

items_comp_scope_exactly_match <- subset(items_comp_scope_exactly, question=="match")

m.scope.exactly.match.0 <- glmer(response_num~modal+(1|id),
                                 glmerControl(optimizer = c("bobyqa")),
                                 family=binomial, data=items_comp_scope_exactly_match)
summary(m.scope.exactly.match.0)


items_comp_scope_less <- subset(items_comp_scope, type=="less")

m.scope.less.0 <- glmer(response_num~question*modal
                           #+ (0 + modal | id) + (0 + modal | item)
                           #+(1|item)
                            +(1|id),
                           family=binomial, data=items_comp_scope_less)
summary(m.scope.less.0)



# #posthoc tests correctd alpha .125
items_comp_scope_modal <- subset(items_comp_scope, modal=="present")

m.scope.modal <- glmer(response_num~type*question
                       #+(1|item)
                       +(1|id),
                       family=binomial, data=items_comp_scope_modal)
summary(m.scope.modal)

items_comp_scope_nomodal <- subset(items_comp_scope, modal=="absent")

m.scope.nomodal <- glmer(response_num~type*question
                         #+(1|item)
                         +(1|id),
                         family=binomial, data=items_comp_scope_nomodal)
summary(m.scope.nomodal)

##Bayes
bm.scope.0.0 <- bglmer(response_num~question*modal*type
                       #+(1|item)
                       +(1|id),
                       family=binomial,
                       data=items_comp_scope, 
                       fixef.prior = normal(cov = diag(9,8)),
                       glmerControl(optimizer = c("bobyqa")))

summary(bm.scope.0.0)

bm.scope.exactly.0 <- bglmer(response_num~question*modal
                          +(1|id),
                          family=binomial,
                          data=items_comp_scope_exactly,
                          #fixef.prior = normal(cov = diag(9,4)),
                          fixef.prior = t(df=7, scale=3),
                          glmerControl(optimizer = c("bobyqa")))
summary(bm.scope.exactly.0)

bm.scope.exactly_mm.0 <- bglmer(response_num~modal
                             +(1|id),
                             family=binomial,
                             data=items_comp_scope_exactly_match,
                             fixef.prior = normal(cov = diag(9,2)),
                             glmerControl(optimizer = c("bobyqa")))
summary(bm.scope.exactly_mm.0)

bm.scope.less.0 <- bglmer(response_num~question*modal
                        +(1|id),
                        family=binomial,
                        data=items_comp_scope_less,
                        fixef.prior = normal(cov = diag(9,4)),
                        glmerControl(optimizer = c("bobyqa")))
summary(bm.scope.less.0)

##save models
save(m.scope.0.0,m.scope.0.1,m.scope.0.3,
     m.scope.0.4,m.scope.0.5,m.scope.0.6,
     m.scope.0.7,m.scope.0.8,m.scope.0.9,
     m.scope.0.all,m.scope.1.8,
     m.scope.exactly.0,m.scope.0.all,
     file = "models.RData") #TODO: save Bayesian models
