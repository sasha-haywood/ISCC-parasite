library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pscl)

alldata <- read_csv("alldata.csv", 
                    col_types = cols(species = col_factor(levels = c("fantail", 
                    "orangethroat", "rainbow")),
                    month = col_factor(levels = c("March", "April",
                    "June", "July", "August", "October")),
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

ggplot(alldata, aes(gill, fill = species)) +
  geom_histogram()

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
final.with = glm(gill ~ alleles + mass + month + species, 
                 alldata, family = "poisson")


# all the orangethroat species have gill=0.  We can build a model with or
# without them and just consider that species separately.

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

# recommended model
final.without = glm(gill ~ alleles + mass + month + species, 
                    data.without.orange, family = "poisson")

# Can compare these possible models
summary(poiss.gill)
summary(final.with)
summary(poiss.gill.without)
summary(final.without)

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

ggplot(alldata, aes(acantho, fill = month)) +
  geom_histogram()
# looks like only species and mass matters

drop1(poiss.acantho, test = "Chi")
# can remove sex
m1 = glm(acantho ~ Gen_freq_all + Gen_freq_species + alleles + 
           mass + month + species, alldata, family = "poisson")
step(m1)
# preferred model does not include alleles but that is what we are testing for,
# so leave it in anyway.

acantho.final = glm(acantho ~ alleles + mass + month + species, alldata, 
                    family = "poisson")
summary(poiss.acantho)
summary(acantho.final)
# alleles is not significant in either model

# These models show some heteroskedasticity.  We can try quasi-poisson or poisson 
# with dispersion. We will try to have a solution to heteroskedasticity by Mon.
