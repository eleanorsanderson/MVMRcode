

#calculating heterogeneity Q statistic for MVMR


##Variables used
#mvmrdata is a dataframe containing a list of snps and the beta and se for each exposure and the outcome 
# here the variables are labeled as
# x1.beta, x1.se <- beta and se for exposure 1
# x2.beta, x2.se <- beta and se for exposure 2
# out.beta, out.se <- beta and se for the outcome

#estbeta1 estimated beta for x1 from the multivariable MR estimation
#estbeta2 estimated beta for x2 from the multivariable MR estimation

#L is the total number of snps in the analysis

seQ <- (mvmrdata$out.se)^2 + (estbeta1^2)*(mvmrdata$x1.se^2) + (estbeta2^2)*(mvmrdata$x2.se^2) 
Q_snp <- ((1/seQ)*((mvmrdata$out.beta-(estbeta1*mvmrdata$x1.beta+estbeta2*mvmrdata$x2.beta))^2))
Q <- sum(Q_snp)
critical_value <- qchisq(0.05,L-2,lower.tail = FALSE)
p_value <- pchisq(Q,L-2,lower.tail = FALSE)

