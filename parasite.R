library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)


alldata <- read_csv("alldata.csv", 
                    col_types = cols(month = col_factor(levels = 
                    c("January", 
                    "February", "March", "April", "May", 
                    "June", "July", "August", "Septemnber", 
                    "October", "November", "December"))))
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
# I think this indiates that we need a hurdle or zero inflated poisson model.

ggplot(alldata, aes(month, gill)) +
  geom_jitter()
ggplot(alldata, aes(month, acantho)) +
  geom_jitter()
# I don't think we have enough months to treat them as anything other than a factor.

# I am not sure we need a mixed model on this.  I don't see what would be the random effect.  
# I also don't think we need to include year.  But I could be wrong.