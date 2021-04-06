library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pscl)
library(DHARMa)
library(lme4)
library(glmmTMB)

alldata <- read_csv("alldata.csv", 
                    col_types = cols(species = col_factor(levels = c("fantail", 
                    "orangethroat", "rainbow")),
                    month = col_factor(levels = c("March", "April",
                    "June", "July", "August", "October", "November")),
                    sex = col_factor(levels = c("F", "M"))))
alldata = alldata[,c(21, 22, 10, 11, 8, 6, 7, 5, 3, 4, 2)]

# Go ahead and remove any observation where we don't have either of the possible
# response variables.
alldata = alldata[!is.na(alldata$gill) | !is.na(alldata$acantho),]


ggplot(alldata, aes(alleles, gill)) +
  geom_jitter()
ggplot(alldata, aes(alleles, acantho)) +
  geom_jitter()
ggplot(alldata, aes(gill)) +
  geom_histogram()
ggplot(alldata, aes(acantho)) +
  geom_histogram()
# Indicates poisson model.

ggplot(alldata, aes(month, gill)) +
  geom_jitter()
ggplot(alldata, aes(month, acantho)) +
  geom_jitter()
# Treat months as a factor.
# Can not treat as time series because we are missing most months.

pairs(alldata[-c(1, 2)])
pairs(alldata[c(6:7)])
# Obviously, mass and length are highly correlated so we have to choose 1.  
# Mass only has 2 missing values and SL has 44, so it makes sense to keep mass.

ggplot(alldata, aes(acantho, fill = sex)) +
  geom_histogram()
ggplot(alldata, aes(gill, fill = sex)) +
  geom_histogram()
# sex probably isn't significant and has lots of missing values.  We will leave it in for now.



##### Poisson model with gill as the response

poiss.gill = glm(gill ~ Gen_freq_all + Gen_freq_species + alleles + 
                   mass +
                   sex + month + species, alldata, family = "poisson")
summary(poiss.gill)

# We can leave all factors in, but if we do model selection, the model 
# improves if we remove sex. Might be good since so many sex values are NA.

drop1(poiss.gill, test = "Chi")

m1 = glm(gill ~ Gen_freq_all + Gen_freq_species + alleles + 
           mass + month + species, alldata, family = "poisson")

# More model selection suggests we can remove Gen_freq_all and Gen_freq_species
step(m1)

# recommended model 
small.with = glm(gill ~ alleles + mass + month + species, 
                 alldata, family = "poisson")



##### all the orangethroat species have gill=0.  We can build a model with or
################# without them and just consider that species separately.

data.without.orange = alldata %>%
  filter(species != "orangethroat")

poiss.gill.without = glm(gill ~ Gen_freq_all + Gen_freq_species + alleles + 
                    mass +
                    sex + month + species, data.without.orange, family = "poisson")
summary(poiss.gill.without)
drop1(poiss.gill.without, test = "Chi")

# Sex can be dropped (lower AIC if it is removed) which is good because lots of 
# missing values for sex.

m1 = glm(gill ~ Gen_freq_all + Gen_freq_species + alleles + 
           mass + month + species, data.without.orange, family = "poisson")
step(m1)

small.without = glm(gill ~ alleles + mass + month + species, 
                    data.without.orange, family = "poisson")


# Can compare these possible models
summary(poiss.gill)
summary(small.with)
summary(poiss.gill.without)
summary(small.without)

# Significant: Alleles, mass, month, species.  
# Not so significant: Gen_freq_all, Gen_freq_species, sex
# It's hard to tell about sex because there are so many missing values, but I think 
# it is safe to remove sex from the model.  If we remove the NAs, sex is only
# significant at .07 level.  Client should decide.  If the only purpose it to show
# that alleles is significant when controlling for everything else, then she might
# want to leave everything in. Alleles is more significant when controlling for 
# everything.




############# Poisson model with acantho as response

ggplot(alldata, aes(acantho, fill = month)) +
  geom_histogram()
# All June and August observations are acantho = 0

# here's what it looks like as a poisson model
poiss.acantho = glm(acantho ~ Gen_freq_all + Gen_freq_species + alleles + mass +
                      sex + month + species, alldata, family = "poisson")
summary(poiss.acantho)
# looks like only species and mass matters

drop1(poiss.acantho, test = "Chi")
# can remove sex
m1 = glm(acantho ~ Gen_freq_all + Gen_freq_species + alleles + 
           mass + month + species, alldata, family = "poisson")
step(m1)
# preferred model does not include alleles but that is what we are testing for,
# so leave it in anyway.

acantho.small = glm(acantho ~ alleles + mass + month + species, alldata, 
                    family = "poisson")
summary(poiss.acantho)
summary(acantho.small)
# alleles is not significant in either model

################# check for homoskedasticity with DHARMa package

plot(simulateResiduals(poiss.gill))
plot(simulateResiduals(poiss.acantho))
# neither are good.  Problems with dispersion, outliers, and heteroskedasticity


### Try with zero inflated models instead

z.gill = glmmTMB(gill ~ Gen_freq_all + Gen_freq_species + alleles + 
                   mass +
                   sex + month + species, alldata,
                   ziformula = ~ 1 , family = "poisson")
plot(simulateResiduals(z.gill))
## This test shows no problems with dispersion or heteroskedasticity or outliers.

#### This is the exact same model, built a package that is more common.  
#### The only benefit of the glmmTMB model is that it can be tested 
#### for heteroskedasticity and other problems with the DHARMa package
z2.gill = zeroinfl(gill ~ Gen_freq_all + Gen_freq_species + alleles + 
                   mass +
                   sex + month + species | 1, alldata,
                   dist = "poisson")


## Looks like we want zero inflated for the gill model after all!

summary(z.gill)

##### zero inflated model with acantho as response

z.acantho = glmmTMB(acantho  ~ Gen_freq_all + Gen_freq_species + alleles + 
                     mass +
                     sex + month + species, alldata,
                   ziformula = ~ 1 , family = "poisson")


plot(simulateResiduals(z.acantho))
summary(z.acantho)
## Looks like we want zero inflated for both.

## The zero-inflated gill model starts by telling us that the odds of getting a
## 0 right off the bat is e^-0.09866:0 = .9061:0
## which could also be stated as a probability of .9061/1.9061 = .48
## The other 52% of the data fits a poisson model where the coefficient for alleles 
## is -.20940 so when alleles increases by 1 unit, the number of gill 
## parasites changes by a factor of e^-20340 = .8111
## Or it decreases by 18.89%

## The zero-inflated acantho model starts by telling us that the odds of getting
## a 0 is e^-0.2750 = .7596 or a probability of .7596/1.7596 = .43
## The rest of the data fits a poisson model where the coefficient for alleles
## is -6.895e-02 so when alleles increases by 1 unit, the number of
## acantho parasites changes by a factor of e^-0.06895 = .9334
## which is not statistically significant.
## Mass and species is much more important.

## If you end up removing sex from the model, you will be able to include more 
## observations and the coefficients will change, but you should be able to 
## calculate the new interpretations.









