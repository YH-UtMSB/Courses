library(MatchIt)

fit1 = lm(
  re78 ~ treat + age + educ + black + hispan + married + nodegree + 
    age*treat + educ*treat + black*treat + hispan*treat + married*treat + nodegree*treat, data = lalonde
)

fit2 = lm(
  re78 ~ treat, data = lalonde
)

summary(fit1)
summary(fit2)