
###################################
#######  ECON 173 Project  ########
# Ben Edelhertz & Marli Stellhorn #
###################################

# Clear the working space
rm(list = ls())

# Set directory
setwd("C:/Users/bened/dev/r")

# Load the packages (all must have been installed)
library(plyr)
library(magrittr)
library(haven)
library(tidyverse)
library(stargazer)
library(plm)
library(ggplot2)

# clustered SEs, clustered on "group"
  
clse = function(reg) { 
    G = length(unique(index(reg,"id")))
    N = length(index(reg,"id"))
    dfa = (G/(G - 1))
    rob = sqrt(diag(dfa*vcovHC(reg, method="arellano", type = "HC1", 
                               cluster = "group")))
    return(rob)
}

##############

### read the data 
op <- as.data.frame(read.csv("opdeath.csv", header=TRUE, sep=",") ) %>% 
    mutate(statedum = as.factor(state),           # state as factor variable
           yeardum = as.factor(year),             # year as factor variable
           time = year-1999)

### Descriptive Statistics
stargazer(op[c("death_rate")], 
          type="text", digits=2, median=TRUE, title="Opioid Data")


########################################
### Medical + Recreational Marijuana ###
########################################

op$med = ifelse((op$state == 'Alaska') | 
                  (op$state == 'Arizona' & op$year >= 2011) |
                  (op$state == 'Arkansas' & op$year >= 2017) |
                  (op$state == 'California') |
                  (op$state == 'Colorado' & op$year >= 2001) |
                  (op$state == 'Connecticut' & op$year >= 2013) |
                  (op$state == 'Delaware' & op$year >= 2012) |
                  (op$state == 'Florida' & op$year >= 2018) |
                  (op$state == 'Hawaii' & op$year >= 2001) |
                  (op$state == 'Illinois' & op$year >= 2014) |
                  (op$state == 'Louisiana' & op$year >= 2016) |
                  (op$state == 'Maine') |
                  (op$state == 'Maryland' & op$year >= 2015) |
                  (op$state == 'Massachusetts' & op$year >= 2013) |
                  (op$state == 'Michigan' & op$year >= 2009) |
                  (op$state == 'Minnesota' & op$year >= 2015) |
                  (op$state == 'Montana' & op$year >= 2005) |
                  (op$state == 'Nevada' & op$year >= 2001) |
                  (op$state == 'New Hampshire' & op$year >= 2014) |
                  (op$state == 'New Jersey' & op$year >= 2011) |
                  (op$state == 'New Mexico' & op$year >= 2008) |
                  (op$state == 'New York' & op$year >= 2015) |
                  (op$state == 'North Dakota' & op$year >= 2017) |
                  (op$state == 'Ohio' & op$year >= 2017) |
                  (op$state == 'Oregon') |
                  (op$state == 'Pennsylvania' & op$year >= 2017) |
                  (op$state == 'Rhode Island' & op$year >= 2007) |
                  (op$state == 'Vermont' & op$year >= 2005) |
                  (op$state == 'Washington') |
                  (op$state == 'Washington DC' & op$year >= 2009) |
                  (op$state == 'West Virginia' & op$year >= 2018) , 1, 0
)


op$rec = ifelse((op$state == 'Alaska' & op$year >= 2015) |
                  (op$state == 'Colorado' & op$year >= 2013) |
                  (op$state == 'Nevada' & op$year >= 2017) |
                  (op$state == 'Oregon' & op$year >= 2015) |
                  (op$state == 'Washington' & op$year >= 2013) , 1, 0
)


###############
### Medical ###
###############

op_med = subset(op, state != 'Alaska' & state != 'Colorado' & state != 'Oregon' & state != 'Nevada' & state != 'Washington')

### Pooled OLS
reg1 <- plm(death_rate ~ med,
             data = op_med, 
             index = c("state","year"), model="pooling")

### State FEs
reg2 = plm(death_rate ~ med,
             data = op_med, 
             index = c("state","year"), 
             model="within", effect="individual")

### Year FEs
reg3 = plm(death_rate ~ med,
             data = op_med, 
             index = c("state","year"), model="fd")
  
### State + Year FEs
reg4 = plm(death_rate ~ med,
             data = op_med, 
             index = c("state","year"), 
             model="within", effect="twoways" )

### State-specific time trend
reg5 = plm(death_rate ~ med + time:statedum,
             data = op_med, 
             index = c("state","year"), 
             model="within", effect="twoways" )

stargazer(reg1, reg2, reg3, reg4, reg5,
            se=list(clse(reg1),clse(reg2),clse(reg3),clse(reg4),clse(reg5)), 
            title="Panel Regressions, Clustered SEs", type="text", 
            keep = c("med","rec"), 
            column.labels=c("OLS", "State-FE", "FD", "St-Yr-FE", "St-Trend"), 
            covariate.labels = c("Medical"),
            dep.var.caption  = "Opioid Deaths per 100,000",
            dep.var.labels   = "Age-Adjusted Rate",
            omit.stat = c("rsq", "f"),
            df=FALSE, digits=2)


####################
### Recreational ###
####################

### Pooled OLS
reg6 <- plm(death_rate ~ rec,
            data = op, 
            index = c("state","year"), model="pooling")

### State FEs
reg7 = plm(death_rate ~ rec,
           data = op, 
           index = c("state","year"), 
           model="within", effect="individual")

### Year FEs
reg8 = plm(death_rate ~ rec,
           data = op, 
           index = c("state","year"), model="fd")

### State + Year FEs
reg9 = plm(death_rate ~ rec,
           data = op, 
           index = c("state","year"), 
           model="within", effect="twoways" )

### State-specific time trend
reg10 = plm(death_rate ~ rec + time:statedum,
           data = op, 
           index = c("state","year"), 
           model="within", effect="twoways" )

stargazer(reg6, reg7, reg8, reg9, reg10,
          se=list(clse(reg6),clse(reg7),clse(reg8),clse(reg9),clse(reg10)), 
          title="Panel Regressions, Clustered SEs", type="text", 
          keep = c("med","rec"), 
          column.labels=c("OLS", "State-FE", "FD", "St-Yr-FE", "St-Trend"), 
          covariate.labels = c("Recreational"),
          dep.var.caption  = "Opioid Deaths per 100,000",
          dep.var.labels   = "Age-Adjusted Rate",
          omit.stat = c("rsq", "f"),
          df=FALSE, digits=2)


#############
### Graph ###
#############


ggplot(subset(op, state == 'Colorado' | state == 'Virginia'), aes(x = year, y = death_rate)) + 
  geom_line(aes(color = state)) +
  scale_color_manual(values = c("chartreuse4", "orange2")) +
  theme_classic() +
  xlab("Year") +
  ylab("Opioid Deaths per 100,000") +
  geom_vline(xintercept = 2013, linetype = 'dashed') +
  geom_vline(xintercept = 2014, linetype = 'dotted')
  


####################################################

####################################################


