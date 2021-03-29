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
alldata = alldata[!is.na(alldata$gill) | !is.na(alldata$acantho),]


ggplot(alldata, aes(alleles, gill)) +
  geom_jitter()
ggplot(alldata, aes(alleles, acantho)) +
  geom_jitter()
ggplot(alldata, aes(gill)) +
  geom_histogram()
ggplot(alldata, aes(acantho)) +
  geom_histogram()
# I think this indicates that we need a hurdle or zero inflated poisson model.

ggplot(alldata, aes(month, gill)) +
  geom_jitter()
ggplot(alldata, aes(month, acantho)) +
  geom_jitter()
# I don't think we have enough months to treat them as anything other than a factor.

# I am not sure we need a mixed model on this.  I don't see what would be the 
# random effect.  
# I also don't think we need to include year.  But I could be wrong.

pairs(alldata[-c(1, 2)])
pairs(alldata[c(6:7)])
# obviously, mass and length are highly correlated.  Mass only has 2 NAs and SL has 
# 44, so it makes sense to keep mass.

z.gill = zeroinfl(gill ~ Gen_freq_all + Gen_freq_species + alleles + mass +
                    sex + month + species, alldata, dist = "poisson")
summary(z.gill)

ggplot(alldata, aes(gill, fill = species)) +
  geom_histogram()

# all the orangethroat species have gill=0.  If we remove those samples and look at 
# them separately, the rest of  the data does not require zero inflated model. 
# Plain poisson works.  Otherwise, species determines the binomial part of the model,
# and then everything except Gen_freq_species is significant in the poisson part.

data.without.orange = alldata %>%
  filter(species != "orangethroat")

poiss.gill.without = glm(gill ~ Gen_freq_all + Gen_freq_species + alleles + 
                    mass +
                    sex + month + species, data.without.orange, family = "poisson")
summary(poiss.gill.without)
drop1(poiss.gill.without, test = "Chi")
# sex can be dropped (lower AIC if it is removed) which is good because lots of 
# missing values for sex

m1 = glm(gill ~ Gen_freq_all + Gen_freq_species + alleles + 
           mass + month + species, data.without.orange, family = "poisson")
step(m1)

# recommended model
final.without = glm(gill ~ alleles + mass + month + species, 
                    data.without.orange, family = "poisson")

# actually, poisson works WITH the orangethroat species as part of it.
poiss.gill = glm(gill ~ Gen_freq_all + Gen_freq_species + alleles + 
                   mass +
                   sex + month + species, alldata, family = "poisson")
summary(poiss.gill)
drop1(poiss.gill, test = "Chi")
# can remove sex
m1 = glm(gill ~ Gen_freq_all + Gen_freq_species + alleles + 
           mass + month + species, alldata, family = "poisson")
step(m1)

# recommended model 
final.with = glm(gill ~ alleles + mass + month + species, 
                 alldata, family = "poisson")


summary(poiss.gill)
summary(final.with)
summary(poiss.gill.without)
summary(final.without)
# significant: Alleles, mass, month, species.  
# Not so significant: Gen_freq_all, Gen_freq_species, sex
# it's hard to tell about sex because there are so many missing values, but I think 
# it is safe to remove sex from the model.  If we remove the NAs, sex is only
# significant at .07 level.  Client should decide.  If the only purpose it to show
# that alleles is significant when controlling for everything else, then she might
# want to leave everything in. Alleles is more significant when controlling for 
# everything.

z.acantho = zeroinfl(acantho ~ Gen_freq_all + Gen_freq_species + alleles + mass +
                    sex + month + species, alldata, dist = "poisson")
summary(z.acantho)

# It looks like month determines the binomial part (and maybe mass), and nothing 
# else is significant. All June and August observations are acantho = 0

ggplot(alldata, aes(acantho, fill = month)) +
  geom_histogram()


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
