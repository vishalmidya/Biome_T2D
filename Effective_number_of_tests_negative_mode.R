# adapted from Li et al. 2012 and Li, 2005 
# https://www.sciencedirect.com/science/article/pii/S0160412021003081?via%3Dihub

ght <- cor(d_logis_dummy[,paste0("neg",seq(1:284))])
et <- eigen(ght)
1/ (sum((et$values>1 + 0)* (et$values - 1)))
# 0.007368629