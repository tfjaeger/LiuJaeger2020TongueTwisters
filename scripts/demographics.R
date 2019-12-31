library(here)
library(readr)
library(dplyr)
library(tidyr)

anonyhits <- read_csv(here('data/hlp_hits_anonymized.csv'))
anonyhits$hitid <- as.factor(anonyhits$hitid)
anonyhits$hittypeid <- as.factor(anonyhits$hittypeid)
anonyhits$assignmentid <- as.factor(anonyhits$assignmentid)
anonyhits$protocol <- as.factor(anonyhits$protocol)
anonyhits$Age <- as.integer(anonyhits$Age)
anonyhits$Sex <- as.factor(anonyhits$Sex)
anonyhits$Race <- as.factor(anonyhits$Race)
anonyhits$Ethnicity <- as.factor(anonyhits$Ethnicity)
anonyhits$Experimenter <- as.factor(anonyhits$Experimenter)
anonyhits$title <- as.factor(anonyhits$title)
anonyhits$UUID <- as.factor(anonyhits$UUID)
anonyhits$Year <- as.factor(anonyhits$Year)
anonyhits$Weekday <- factor(anonyhits$Weekday, levels = c("Sun", "Sat", "Fri", "Thu", "Wed", "Tue", "Mon"), ordered = TRUE)
anonyhits$Hour <- as.factor(anonyhits$Hour)
anonyhits$Experiment <- as.factor(anonyhits$Experiment)

# unique(anonyhits$Experiment)
# there was another "Tongue Twister" experiment in 2016. Excluding that here.
TT <- anonyhits %>% 
  filter(Experiment %in% c("TT")) %>%
  select(UUID, Race, Sex, Ethnicity, Age, Year) %>%
  mutate(Sex = Sex %>% replace_na("Unknown"),
         Race = Race %>% replace_na("Unknown or Not Reported"),
         Ethnicity = Ethnicity %>% replace_na("Unknown")) %>%
  droplevels()
summary(TT)
length(unique(TT$UUID))

round(prop.table(table(TT$Ethnicity))*100,1)
round(prop.table(table(TT$Race))*100,1)
round(prop.table(table(TT$Sex))*100,1)

summary(TT$Age)
hist (TT$Age, breaks = 100)
sd(TT$Age, na.rm = T)

