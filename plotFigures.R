# Calculate MSY, Fmsy based on life histories and selectivity from RAM 

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

# Load in the relevant species 
marine.df <- read.csv('marinecatches.csv')

feed.history <- read.csv('Feed projection.csv') # Feed projections 
#feed.projection <- read.csv('ALL_ForageFish_results.csv')
feed.projection <- read.csv('New_Estimates2.csv')




# Sum Catches over year and size 
marine.df.sum <- marine.df %>%
  group_by(Year, size) %>%
  summarise(Catch = sum(Catch, na.rm = T))

sumCatch <- marine.df.sum %>%
  group_by(Year) %>%
  summarise(Catch = sum(Catch, na.rm = T))


p1 <- ggplot(marine.df.sum, aes(x= Year, y = Catch*1e-6, group= size, color = size))+geom_line()+theme_classic()+geom_line(data=)+
  scale_y_continuous('Catch (million tonnes)')
p1
## Plot the fraction of forage fish 
p2 <- ggplot(data= marine.df.sum[marine.df.sum$size == 'Small',], aes(x = Year, y= Catch/sumCatch$Catch, group =1) )+geom_line()+theme_classic()+
  scale_y_continuous('Fraction forage fish of total catch')
p2

df.Forage <- marine.df[marine.df$size == 'Small',] # Subset of catches that are forage fish

df.Forage$stock <- paste(df.Forage$Species,df.Forage$Fishingarea, sep = '-') # Unique FAO stock

sp <- unique(df.Forage$Species)

# Order by total catch in year 2012 
df.2010 <- df.Forage[df.Forage$Year == 2012,]
# Find the unique species (or genus)
df.2010 <- df.2010 %>%
  group_by(Species) %>%
  summarise(Catch = sum(Catch, na.rm = T))

df.2010 <- df.2010[order(df.2010$Catch, decreasing = T),]

# Calculate the forage fish that contribute the most to total catch
df.2010$frac <- df.2010$Catch/sum(df.2010$Catch, na.rm = T)
df.2010$csum <- cumsum(df.2010$frac)

# Create table for export 
df.2010exp <- data.frame(Species = df.2010$Species[1:20], Catch = df.2010$Catch[1:20])

#write.table(df.2010exp, file = 'forage20exp.csv', sep = '\t',row.names = F)

sumForage <- df.Forage %>%
  group_by(Year) %>%
  summarise(Catch = sum(Catch, na.rm =T))


Forage <- df.Forage %>%
  group_by(Year,Species) %>%
  summarise(Catch = sum(Catch, na.rm = T), wmax = mean(wmax, na.rm = T))

Forage <- Forage[Forage$Year > 1961,]

Forage <- Forage[order(Forage$Catch, decreasing = T),]


nspecies <- length(unique(Forage$Species))
species <- unique(Forage$Species)
Mcatch <- matrix(NA, nspecies)
SDcatch <- matrix(NA, nspecies)
SEcatch <- matrix(NA, nspecies)

for (i in 1:nspecies){
  idx <- which(Forage$Species[Forage$Year > 1980] == species[i])
  Mcatch[i] <- mean(Forage$Catch[idx], na.rm =T)
  SDcatch[i] <- sd(Forage$Catch[idx], na.rm =T)
  SEcatch[i] <- SDcatch[i]/sqrt(length(idx))
  
}
# Plot from 85 onwards, mean after 

year <- seq(2012,2050)

# Landings increase over the last 60 years 
newdata <- data.frame(year = year, Cmedian = mean(sumForage$Catch[sumForage$Year > 1980], na.rm = T),
                      SE = sd(sumForage$Catch))

# Plot in ggplot 
sumForage$scenario <- 'Historic Catch'
sumForage$SE <- NA

newdata <- data.frame(Year = year, Catch = mean(sumForage$Catch[sumForage$Year > 1980], na.rm = T),
                      SE = NA, scenario = 'avg. Historic catch > 1980')
newdata$SE <- sd(sumForage$Catch[sumForage$Year > 1980], na.rm = T)/sqrt(nspecies)
newdata$SD <- sd(sumForage$Catch[sumForage$Year > 1980], na.rm = T)

# Add the projection data 
tmp <- rep(NA, length(feed.projection$Mean)*2)
tmp[seq(2,length(feed.projection$Mean)*2, by = 2)] <- feed.projection$Mean
tmp[seq(1,length(feed.projection$Mean)*2-1, by = 2)] <- feed.projection$Mean[1]
proj <- data.frame(Year = rep(c(2012,2050),length(feed.projection$Mean)*2), 
                   Catch = tmp, 
                   SE = rep(feed.projection$SE,each = 2), 
                   scenario = rep(feed.projection$Scenario,each = 2))

## Remove some of the scenarions keep for later
# scns <- unique(proj$scenario)
# proj <- proj[-which(proj$scenario == scns[5]),]
# proj <- proj[-which(proj$scenario == scns[6]),]
# proj <- proj[-which(proj$scenario == scns[1]),]

# Make exponential growth towards 2050
scns <- unique(proj$scenario)
yr <- 2012:2050

df.plot <- data.frame(
  scenario = rep(scns, each = length(yr)), 
  Year = rep(yr, length(scns)), 
  Catch = NA, 
  SE = NA)

for (i in 1:length(scns)){
  df.tmp <- proj[proj$scenario == scns[i],]
  
  y.ml <- lm(log(Catch) ~ Year, data = df.tmp)
  y.pred <- exp(predict(y.ml, newdata = data.frame(Year = yr)))
  
  df.plot[df.plot$scenario == scns[i],]$Catch <- y.pred
  df.plot[df.plot$scenario == scns[i],]$SE <- df.tmp$SE[1]
}

pelfish <- data.frame(Year = feed.history$Year,
                      Catch = feed.history$Pelagic.Fish, scenario = 'Feed use', SE = NA)


# sumForage.plot <- rbind(sumForage,newdata, proj, pelfish) # for plotting
# sumForage.feed <- rbind(proj,pelfish)

df.msy.all.2 <- read.csv('CatchMSY_Nis_FAOAreas_Rtest.csv', sep =',')
df.msy.all.2$stock <- paste(df.msy.all.2$Species,df.msy.all.2$Area)
df.msy.all.2 <- df.msy.all.2[-which(duplicated(df.msy.all.2$stock) == 1),]# remove dupes  


F50 <- read.csv('F_50Fmsy_median.csv')
# sum of MSY 
SE <- c(sum(df.msy.all$quantile05msy),sum(df.msy.all$quantilemsyr95))

df.tmp = data.frame(Catch = c(sumForage$Catch[length(sumForage$Catch)], sum(df.msy.all$median.msy.)),
                    Year = c(2012,2050))

y.ml <- lm(log(Catch) ~ Year, data = df.tmp)
y.pred <- exp(predict(y.ml, newdata = data.frame(Year = yr)))
y.pred <- rep(df.tmp$Catch[2], length(yr))

# Find the difference between median and SE's
SE.diff <- c(sum(df.msy.all$median.msy.)-SE[1],SE[2]-sum(df.msy.all$median.msy.)) 

MSY <- data.frame(Catch = y.pred,
                  Year = yr,scenario = 'MSY', SE.min = y.pred-SE.diff[1], SE.max = y.pred+SE.diff[2])


# Plot the feed scenarios and the catch scenarios 
# Find the order of scenarios 
scen.order <- proj[proj$Year == 2050,]
scen.order <- scen.order[order(scen.order$Catch, decreasing = T),]
brks <- scen.order$scenario
##
pelfishMean <- data.frame(Year = year, 
                          Catch = rep(mean(pelfish$Catch[pelfish$Year > 1980]),length(year)),
                          SD = rep(sd(pelfish$Catch[pelfish$Year > 1980]),length(year)),
                          scenario = 'Historical feed')


p1 <- ggplot(df.plot, aes(x = Year, y = Catch*1e-6, color = scenario))+
  geom_line(linetype = 'dashed')+theme_classic()+ 
  geom_ribbon(aes(ymin = (Catch-SE)*1e-6, ymax = (Catch+SE)*1e-6, fill = scenario), alpha = 0.15, linetype = 0)+
  geom_line(data = pelfish, col = alpha('red3', alpha = 0.8), linetype = 'dashed')+
  #  geom_line(data = pelfishMean, col = alpha('red3', alpha = 0.8), linetype = 'dashed')+
  #geom_ribbon(data = pelfishMean, fill = alpha('red3'), alpha = 0.3, linetype = 0, 
  #            aes(ymin = (Catch-SD)*1e-6, ymax = (Catch+SD)*1e-6 ))+
  geom_line(data= sumForage, col = alpha('dodgerblue', alpha = 1),aes(y = Catch*1e-6*.9))+
  geom_line(data = newdata, col = alpha('dodgerblue', alpha = 1), aes(y = Catch*1e-6*0.9))  + 
  # geom_ribbon(data = newdata, aes(x = Year, 
  #                                 ymin = (Catch-SD)*1e-6*0.9, 
  #                                 ymax = (Catch+SD)*1e-6*0.9), 
  #                                 fill = 'dodgerblue', alpha= 0.3, linetype = 0)+
  geom_vline(xintercept = 2012, linetype = 2, color = 'black')+
  geom_line(data = MSY, col = 'royalblue3', aes(y = Catch*1e-6*0.9))+
  #   geom_ribbon(data = MSY, aes(x = Year, ymin = SE.min*1e-6*0.9, ymax = SE.max*1e-6*0.9), fill = 'royalblue4', alpha= 0.3, linetype = 0)+
  scale_y_continuous('Biomass (million tonnes)')+
  annotate('text',x = 2025, y = 45, label = 'Projected max supply', colour = 'royalblue3', size = 2)+
  annotate('text',x = 2035, y = 12, label = 'Avg. historical supply (>1980)', colour = 'dodgerblue', size = 2)+
  annotate('text',x = 1995, y = 35, label = 'Historical supply', colour = 'dodgerblue', size = 2)+
  annotate('text',x = 1990, y = 10, label = 'Historical feed use', colour = 'red3', size = 2)+
  scale_color_discrete("", breaks = brks)+
  guides(fill = FALSE, color=guide_legend(override.aes=list(fill=NA)))+
  theme(legend.position=c(.18,.8),text = element_text(size=8),
        legend.text = element_text(size=7),legend.direction = "vertical",
        plot.margin=unit(c(0,0.3,0.2,0.2), "cm"),
        legend.background = element_rect(fill="transparent"),
        legend.key.width=unit(0.7,"line"),
        legend.key.height=unit(0.7,"line"))

windows(width = 3.15*2*2/3, height = 3.5)  
p1
#legend.box.background = element_rect()
print(paste('MSY Increase over 2012 catch =', (MSY$Catch[1]/sumForage$Catch[length(sumForage$Catch)]-1)*100), '%')
# Biomass
print(paste('MSY Increase over 2012 catch =', (MSY$Catch[1]-sumForage$Catch[length(sumForage$Catch)])*1e-6, 'million tonnes'))

print(paste('MSY Increase over 2012 catch =', (MSY$Catch[1]/newdata$Catch[1]-1)*100, '%'))
print(paste('MSY Increase over 2012 catch =', (MSY$Catch[1]-newdata$Catch[1])*1e-6, 'million tonnes'))

# Forage fish reductions 

print(paste('50 % decrease = ', sum(df.msy.all$mediank)*0.2*1e-6*0.9)) # 0.9 is availability for consumption
print(paste('20 % decrease = ', MSY$Catch[1]*0.8*1e-6*0.9))

# Calculate potential increase in catch
#  scale_fill_manual("Catch", values = 1:9)
# Difference between use and catch 

cairo_pdf('C:/Users/Nis/Dropbox/UW/Aquaculture Halley/Figures/TotalUSE.pdf', width = 3.15*2*2/3, height = 3.5)
p1
dev.off()

jpeg('C:/Users/Nis/Dropbox/UW/Aquaculture Halley/Figures/TotalUSE.jpg', width = 3.15*2*2/3, height = 3.5, unit = 'in', res = 1200)
p1
dev.off()

# Print the years where they cross mean supply and MSY
scenarios <- unique(df.plot$scenario)

df.yr.cross <- data.frame(scenario = scenarios, yr.msy = NA, yr.msy.c = NA, yr.med = NA, yr.med.c = NA, yr.msy20 = NA, yr.msy20.c = NA,
                          yr.mid.min = NA, yr.mid.max = NA,
                          yr.msy.min = NA, yr.msy.max = NA)

for(i in 1:length(scenarios)){
  msy <- MSY$Catch[1]*.9
  med.c <- newdata$Catch[1]*.9
  med.c.min <- newdata$Catch[1]*0.9-newdata$SD[1]*0.9
  med.c.max <- newdata$Catch[1]*0.9+newdata$SD[1]*0.9
  msy20 <- MSY$Catch[1]*.9*.8 # Leave 20% less than now for the hangray preds
  
  msy.min <- MSY$SE.min*0.9
  msy.max <- MSY$SE.max*0.9
  
  tmp <- df.plot[df.plot$scenario == scenarios[i],]
  
  df.yr.cross$yr.msy[i] <-tmp$Year[which.min((msy - tmp$Catch)^2)]
  df.yr.cross$yr.med[i] <-tmp$Year[which.min((med.c - tmp$Catch)^2)]
  df.yr.cross$yr.msy20[i] <-tmp$Year[which.min((msy20 - tmp$Catch)^2)]
  df.yr.cross$yr.mid.min[i] <-tmp$Year[which.min((med.c.min - tmp$Catch)^2)]
  df.yr.cross$yr.mid.max[i] <-tmp$Year[which.min((med.c.max - tmp$Catch)^2)]
  
  
  df.yr.cross$yr.msy.c[i] <-tmp$Catch[which.min((msy - tmp$Catch)^2)]
  df.yr.cross$yr.med.c[i] <-tmp$Catch[which.min((med.c - tmp$Catch)^2)]
  df.yr.cross$yr.msy20.c[i] <-tmp$Catch[which.min((msy20 - tmp$Catch)^2)]
  
  df.yr.cross$yr.msy.min[i] <- tmp$Year[which.min((msy.min - tmp$Catch)^2)]
  df.yr.cross$yr.msy.max[i] <- tmp$Year[which.min((msy.max - tmp$Catch)^2)]
  
  
}
# df.yr.cross$yr.msy[df.yr.cross$yr.msy == 2050] <- NA
# df.yr.cross$yr.med[df.yr.cross$yr.med == 2050] <- NA
df.yr.cross$yr.msy[df.yr.cross$yr.msy == 2012] <- NA
df.yr.cross$yr.med[df.yr.cross$yr.med == 2012] <- NA
#df.yr.cross$yr.yr.msy20[df.yr.cross$yr.msy20 == 2050] <- NA
# df.yr.cross$yr.yr.mid.max[df.yr.cross$yr.mid.max == 2050] <- NA
# df.yr.cross$yr.mid.min[df.yr.cross$yr.mid.min == 2050] <- NA
df.yr.cross$yr.msy20[df.yr.cross$yr.msy20 == 2012] <- NA
df.yr.cross$yr.med[df.yr.cross$yr.med == 2012] <- NA
df.yr.cross$yr.mid.max[df.yr.cross$yr.mid.max == 2012] <- NA
df.yr.cross$yr.mid.min[df.yr.cross$yr.mid.min == 2012] <- NA

df.yr.cross

jpeg('C:/Users/Nis/Dropbox/UW/Aquaculture Halley/Figures/TotalUSE_yrs.jpg', width = 3.15*2*2/3, height = 3.5, unit = 'in', res = 1200)

p1 + geom_point(data = df.yr.cross, aes(x = yr.msy, y = yr.msy.c*1e-6))+
  geom_point(data = df.yr.cross, aes(x = yr.med, y = yr.med.c*1e-6))

dev.off()
write.table(df.yr.cross, 'omega_deficit_years.csv', row.names = F, sep = ',')

# Add the barplot to fig 2 

df.spec <- read.csv('New_Current_Taxon_Fig.csv')
df.spec <- df.spec[order(df.spec$Current_Mean),]
df.spec$Group <- factor(df.spec$Group, levels = df.spec$Group[order(df.spec$Current_Mean, decreasing = T)])

p2 <- ggplot(data = df.spec, aes(y = Current_Mean_scaled, x= as.factor(Group)))+geom_col(aes(fill = ID))+theme_classic()+coord_flip()+
  scale_y_continuous('Biomass (million tons)')+ scale_x_discrete('')+
  geom_errorbar(
    aes(ymin = Current_Mean_scaled-Current_SD_scaled, 
        ymax = Current_Mean_scaled+Current_SD_scaled), size = 0.01, width = 0.3)+
  theme(legend.position='none',
        plot.margin=unit(c(0.1,0.1,0.2,-1), "cm"),text = element_text(size=8))

grid.arrange(p1,p2, ncol = 2)

# 
cairo_pdf('C:/Users/Nis/Dropbox/UW/Aquaculture Halley/Figures/speccontr.pdf', width = 3.15, height = 3.15/2)
p2
dev.off()

jpeg('C:/Users/Nis/Dropbox/UW/Aquaculture Halley/Figures/speccontr.jpg', width = 3.15*2, height = 3.15, unit = 'in', res = 1200)
p2
dev.off()



### PLot the two figures in the same graph 
cairo_pdf('C:/Users/Nis/Dropbox/UW/Aquaculture Halley/Figures/speccontr.pdf', width = 3.15*2, height = 3.15)
grid.arrange(p1,p2, ncol = 2)
dev.off()
# 
# vp <- viewport(width =0.1, height = 0.1, x = 1950,
#                y = unit(0.7, "lines"), just = c("left",
#                                                 "top"))
# 
# 
# full <- function() {
#   print(p1)
#   print(p2, vp = vp)
#   
# }
# 
# 
# full()
# 

# Difference between use and catch 



yr <- 1961:2012
diff.c <- sumForage$Catch[sumForage$Year > (min(pelfish$Year)-1)]-pelfish$Catch[pelfish$Year < (max(sumForage$Year)+1)] 
plot(yr, diff.c, type = 'l')
ml1 <- lm(diff.c ~ yr) 
abline(ml1)
summary(ml1)


yr <- 1980:2012
diff.c <- pelfish$Catch[pelfish$Year < 2013 & pelfish$Year> 1979]/sumForage$Catch[sumForage$Year > 1979]

cff <- coef(ml1 <- lm(diff.c ~ yr)) 

p3 <- ggplot(data.frame(yr = yr, diff.c = diff.c), aes(x = yr, y = diff.c))+geom_line()+geom_point()+theme_classic()+
  geom_smooth(method='lm', fill = 'blue', alpha = 0.2)+scale_y_continuous('fraction of total catch')+scale_x_continuous('Year')
p3

jpeg(filename = 'fractionofcatch.jpg', width = 16,height = 16, unit = 'cm', res = 600)
p3
dev.off()


# How large of a fraction of the catch is the total msy calculation in 2012 


df.tmp <- df.Forage[df.Forage$Year == 2012,]
# Remove dupes 
df.tmp <- df.tmp[-which(duplicated(df.tmp$stock) == 1),] # Need to find those dupes later 


df.msy.all$stock2 <- paste(df.msy.all$Species,df.msy.all$Area, sep = '-')
nstocks <- length(df.tmp$stock)

ixvec <- matrix(0, nstocks)
for (i in 1:nstocks){
  ixtmp <- which(df.tmp$stock[i] == df.msy.all$stock2)  
  
  if(length(ixtmp) > 0){
    ixvec[i] = 1
  }
}


# Sum UP! 
sum(df.tmp[ixvec == 1,]$Catch, na.rm = T)/sum(df.tmp$Catch, na.rm = T)
# Number of speices 
length(unique(df.msy.all$Species))

