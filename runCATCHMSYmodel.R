# Run catch MSY with FAO data 
# Majority of code from Martell and Froese 2013, Fish and Fisheries
## Load packages 
library(dplyr)
library(ggplot2)
library(rfishbase)
set.seed(1000)

## Read FAO data 
df.fao <- read.csv('marinecatches.csv')
head(df.fao)

idx <- which(df.fao$size == 'Small')
df.forage <- df.fao[idx,]
areas <- unique(df.fao$Fishingarea)

# Plot all the forage catches 
df.tot <- df.forage %>%
  group_by(Year) %>%
  summarise(Catch = sum(Catch, na.rm = T))

head(df.tot)

ggplot(df.tot, aes(x = Year, y = Catch*1e-6))+geom_line()+theme_classic()+scale_y_continuous('Total Catch')

sp <- unique(df.forage$Species)

# Create a new data frame with unique stocks per fishing area 
df.forage$stock <- paste(df.forage$Species,df.forage$Fishingarea, sep = '-')
head(df.forage)

stock <- unique(df.forage$stock)
stock <- rev(stock) # ORdered with smallest catch in the beginning

outfile  <- "CatchMSY_FAOAreas.csv"

for(i in 1:length(stock)) {
  df.temp <- df.forage[df.forage$stock == stock[i],]
  
  df.temp <- df.temp[order(df.temp$Year),]
  species <- df.temp$Species[1]
  
  yr   <- df.temp$Year
  ct   <- df.temp$Catch  
  ct[is.na(ct)] = 0 # NA is assumed zero for simplicity
  
  nyr  <- length(yr)    ## number of years in the time series
  
  if (sum(ct) > 100){ # dont worry about the very small stocks
  
  start_r     <- c(0.05,2)  # Initial r values - we use a big span to include all species
  start_k     <- c(max(ct),100*max(ct)) ## default for upper k e.g. 100 * max catch
  
  startbio    <- if(ct[1]/max(ct, na.rm =T) < 0.5) {c(0.5,0.9)} else {c(0.3,0.6)} ## use for batch processing
  interyr 	<- yr[2]   ## interim year within time series for which biomass estimate is available; set to yr[2] if no estimates are available
  interbio 	<- c(0, 1) ## biomass range for interim year, as fraction of k; set to 0 and 1 if not available
  finalbio    <- if(ct[nyr]/max(ct, na.rm =T) > 0.5) {c(0.3,0.7)} else {c(0.01,0.5)} ## use for batch processing
  n           <- 30000  # Expanding this number takes a lot of time
  
  startbt     <- seq(startbio[1], startbio[2], by = 0.05) ## apply range of start biomass in steps of 0.05	
  parbound <- list(r = start_r, k = start_k, lambda = finalbio)
  
  ## FUNCTIONS
  runModel	<- function(theta)
  {
    with(as.list(theta), {  ## for all combinations of ri & ki
      bt=vector()
      ell = 0  ## initialize ell
      for (j in startbt)
      {
        if(ell == 0) 
        {
          bt[1]=j*k  ## set biomass in first year
          for(i in 1:nyr) ## run the Schaefer model
          {
              bt[i+1]=(bt[i]+r*bt[i]*(1-bt[i]/k)-ct[i]) ## calculate biomass as function of previous year's biomass plus net production minus catch
          }
          
          #Bernoulli likelihood
          ell = 0
          if(bt[nyr+1]/k>=lam1 && bt[nyr+1]/k <=lam2 && min(bt) > 0 && max(bt) <=k && 
             bt[which(yr==interyr)]/k>=interbio[1] && bt[which(yr==interyr)]/k<=interbio[2]) 
            ell = 1
        }	
      }
      return(list(ell=ell))
      
      
    })
  }
  
  sraMSY	<-function(theta, N)
  {
    #This function conducts the stock reduction analysis
    
    with(as.list(theta), 
         {
           ri = exp(runif(N, log(r[1]), log(r[2])))  ## get N values between r[1] and r[2], assign to ri
           ki = exp(runif(N, log(k[1]), log(k[2])))  ## get N values between k[1] and k[2], assing to ki
           itheta=cbind(r=ri,k=ki, lam1=lambda[1],lam2=lambda[2]) ## assign ri, ki, and final biomass range to itheta
           M = apply(itheta,1,runModel) ## call Schaefer function with parameters in itheta
           i=1:N
           ## prototype objective function
           get.ell=function(i) M[[i]]$ell
           ell = sapply(i, get.ell) 
           return(list(r=ri,k=ki, ell=ell))	
         })
  }
  
  ## MAIN
  R1 = sraMSY(parbound, n)  
  
  ## Get statistics on r, k, MSY and determine new bounds for r and k
  r1 	<- R1$r[R1$ell==1]
  k1 	<- R1$k[R1$ell==1]
  msy1  <- r1*k1/4 
  mean_msy1 <- exp(mean(log(msy1))) 
  max_k1a  <- min(k1[r1<(1.1*parbound$r[1])]) ## smallest k1 near initial lower bound of r
  max_k1b  <- max(k1[r1*k1/4<mean_msy1]) ## largest k1 that gives mean MSY
  max_k1 <- if(max_k1a < max_k1b) {max_k1a} else {max_k1b}
  
  if(length(r1)<10) {
    cat("Too few (", length(r1), ") possible r-k combinations, check input parameters","\n")
  }
  
  if(length(r1)>=10) { # Only run with sufficient r values
    
    ## set new upper bound of r to 1.2 max r1
    parbound$r[2] <- 1.2*max(r1)
    ## set new lower bound for k to 0.9 min k1 and upper bound to max_k1 
    parbound$k 	  <- c(0.9 * min(k1), max_k1)
    
    
    cat("First MSY =", format(1000*mean_msy1, digits=3),"\n")
    cat("First r =", format(exp(mean(log(r1))), digits=3),"\n")
    cat("New upper bound for r =", format(parbound$r[2],digits=2),"\n")	
    cat("New range for k =", format(1000*parbound$k[1], digits=3), "-", format(1000*parbound$k[2],digits=3),"\n")
    
    
    ## Repeat analysis with new r-k bounds
    R1 = sraMSY(parbound, n)
    
    ## Get statistics on r, k and msy
    r = R1$r[R1$ell==1]
    k = R1$k[R1$ell==1]
    msy = r * k / 4
    mean_ln_msy = mean(log(msy))
    
    ## plot MSY over catch data - uncomment to plot
    # fig.name <- paste(stock[i], sep = '-')
    # fig.name <- paste(fig.name,".jpg", sep = '')
    # 
    # 
    # # Plot the data to manually check MSY values
    # 
    # jpeg(filename = paste("FiguresMSY_testR",fig.name,sep = "/"), width = 16, height = 12, units = 'cm', res = 100)
    
    par(mfcol=c(2,3))
    plot(yr, ct, type="l", ylim = c(0, max(ct)), xlab = "Year", ylab = "Catch", main = species)
    abline(h=exp(mean(log(msy))),col="red", lwd=2)
    abline(h=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
    abline(h=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
    
    hist(r, freq=F, xlim=c(0, 1.2 * max(r)), main = "")
    abline(v=exp(mean(log(r))),col="red",lwd=2)
    abline(v=exp(mean(log(r))-2*sd(log(r))),col="red")
    abline(v=exp(mean(log(r))+2*sd(log(r))),col="red")
    
    plot(r1, k1, xlim = start_r, ylim = start_k, xlab="r", ylab="k (1000t)")
    
    hist(k, xlim=c(0, 1.2 * max(k)), xlab="k (1000t)", main = "")
    abline(v=exp(mean(log(k))),col="red", lwd=2)	
    abline(v=exp(mean(log(k))-2*sd(log(k))),col="red")
    abline(v=exp(mean(log(k))+2*sd(log(k))),col="red")
    
    
    plot(log(r), log(k),xlab="ln(r)",ylab="ln(k)")
    abline(v=mean(log(r)))
    abline(h=mean(log(k)))
    abline(mean(log(msy))+log(4),-1, col="red",lwd=2)
    abline(mean(log(msy))-2*sd(log(msy))+log(4),-1, col="red")
    abline(mean(log(msy))+2*sd(log(msy))+log(4),-1, col="red")
    
    hist(msy, freq=F, xlim=c(0, 1.2 * max(msy)), xlab="MSY (1000t)",main = "")
    abline(v=exp(mean(log(msy))),col="red", lwd=2)
    abline(v=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
    abline(v=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
    
    #dev.off()

    ## Write results into outfile, in append mode (no header in file, existing files will be continued)
    output = data.frame(paste(species,stock[i], sep = '-'), startbio[1], startbio[2], 
                        interbio[1], interbio[2], finalbio[1], finalbio[2], 
                        min(yr), max(yr), max(ct), ct[1], 
                        ct[nyr], length(r), exp(mean(log(r))), 
                        sd(log(r)), min(r), quantile(r,0.05), 
                        quantile(r,0.25), median(r), quantile(r,0.75), 
                        quantile(r,0.95), max(r), exp(mean(log(k))), sd(log(k)), 
                        min(k), quantile(k, 0.05), quantile(k, 0.25), median(k), 
                        quantile(k, 0.75), quantile(k, 0.95), max(k), 
                        exp(mean(log(msy))), sd(log(msy)), min(msy), 
                        quantile(msy, 0.05), quantile(msy, 0.25), median(msy), 
                        quantile(msy, 0.75), quantile(msy, 0.95), max(msy)) 
    
    write.table(output, file = outfile, append = TRUE, sep = ";", dec = ".", row.names = FALSE, col.names = FALSE)
    
  }
  print(i) # How many have we run? 
  }
} 

