
library(rBBS)
library(RPresence)
library(BBSclim)

### Import raw 10-stop BBS data
bbs <- GetRouteData()

### Import weather data
weather <- GetWeather()

### Import route information
routes <- GetRoutes()

### Get count data for WOTH
alpha <- "lowa"

spp_AOU <- GetAOU(alpha)

GetSppCounts(AOU = spp_AOU)


### Read raw BBS counts
spp_counts <- read.csv(paste("inst/output/spp_counts", paste(alpha, "counts.csv", sep = "_"), sep = "/"))
spp_counts2 <- dplyr::filter(spp_counts, Longitude > -110)
buffer_BBS(spp_count = spp_counts2)


### Load rasters containing bioclim estimates
data("NA_biovars")

### Read buffered BBS counts
spp_buff <- read.csv(paste("inst/output/buffered_counts",
                           paste(alpha, "buff.csv", sep = "_"), sep = "/"))

### Extract annual bioclim estimates for Wood Thrush BBS routes
GetBioVars(counts = spp_buff)


### Save occupancy and climate data as .pao file for input into Presence
spp_clim <- read.csv(paste("inst/output/spp_clim",
                           paste(alpha, "clim.csv", sep = "_"), sep = "/"))

library(RPresence)
write_pao(counts = spp_buff, clim = spp_clim, alpha = alpha)

### Create design matrices

spp_pao <- read.pao(paste("inst/output/pao", paste(alpha, ".pao", sep = ""), sep = "/"))

psi_mods <- GetPsiMods()

RunPsiMods(pao = spp_pao, alpha = alpha)



###	read the output of the file we just ran
list.out <- list.files('inst/output/pres/','.*.out$')
model.id <- grep(modname, list.out)
model.out <- scan(paste0('inst/output/pres/', list.out[model.id]), what='character', sep='\n', quiet=T)

jj <- grep('std.error', model.out)
jj2 <- grep('Individual Site estimates of <psi>', model.out)
betas=model.out[(jj+1):(jj2-1)]				# works for het and plain?
coefs <- as.numeric(substr(betas, 41,50))
std.er <- as.numeric(substr(betas, 54,63))

psi.keep <- th0.keep <- th1.keep <- gam.keep <- eps.keep <- p1.keep <- p2.keep <- stop.keep <- NULL
### drop covariates with negative standard errors
if(num.betas[1]>1) {
  psi.keep 	<- !is.na(std.er[2:num.betas[1]])
  #psi.keep[1]=FALSE; psi.keep[4]=FALSE		# just for testing
  psi.cov	<- drop.bad.coef(psi.keep, psi.cov)
}
if(num.betas[2]>1) {
  th0.keep 	<- !is.na(std.er[(2+num.betas[1]):sum(num.betas[1:2])])
  th0.cov	<- drop.bad.coef(th0.keep, th0.cov)
}
if(num.betas[3]>1) {
  th1.keep 	<- !is.na(std.er[(2+sum(num.betas[1:2])):sum(num.betas[1:3])])
  th1.cov	<- drop.bad.coef(th1.keep, th1.cov)
}
if(num.betas[4]>1) {
  gam.keep 	<- !is.na(std.er[(2+sum(num.betas[1:3])):sum(num.betas[1:4])])
  gam.cov	<- drop.bad.coef(gam.keep, gam.cov)
}
if(num.betas[5]>1) {
  eps.keep 	<- !is.na(std.er[(2+sum(num.betas[1:4])):sum(num.betas[1:5])])
  eps.cov	<- drop.bad.coef(eps.keep, eps.cov)
}

if(num.betas[6]>1) {
  p1.keep   	<- p2.keep  <- !is.na(std.er[(2+sum(num.betas[1:5])):sum(num.betas[1:6])])
  if(het==T)      p2.keep  <- !is.na(std.er[(2+sum(num.betas[1:6])):(sum(num.betas[1:6])+num.betas[6])])
  p1.cov	<- drop.bad.coef(p1.keep & p2.keep, p1.cov)
}



### drop covariates with large coefficients (if no neg SE)
### change to drop cov with large SE (if no neg SE, i.e., keep==1)
if (mean(c(psi.keep, th0.keep, th1.keep, gam.keep, eps.keep, p1.keep, p2.keep, 1))==1) {

  large=2
  if(num.betas[1]>1) {
    psi.keep 	<- abs(std.er[2:num.betas[1]])<large
    #psi.keep[1]=FALSE; psi.keep[4]=FALSE		# just for testing
    psi.cov	<- drop.bad.coef(psi.keep, psi.cov)
  }
  if(num.betas[2]>1) {
    th0.keep 	<- abs(std.er[(2+num.betas[1]):sum(num.betas[1:2])])<large
    th0.cov	<- drop.bad.coef(th0.keep, th0.cov)
  }
  if(num.betas[3]>1) {
    th1.keep 	<- abs(std.er[(2+sum(num.betas[1:2])):sum(num.betas[1:3])])<large
    th1.cov	<- drop.bad.coef(th1.keep, th1.cov)
  }
  if(num.betas[4]>1) {
    gam.keep 	<- abs(std.er[(2+sum(num.betas[1:3])):sum(num.betas[1:4])])<large
    gam.cov	<- drop.bad.coef(gam.keep, gam.cov)
  }
  if(num.betas[5]>1) {
    eps.keep 	<- abs(std.er[(2+sum(num.betas[1:4])):sum(num.betas[1:5])])<large
    eps.cov	<- drop.bad.coef(eps.keep, eps.cov)
  }
  if(num.betas[6]>1) {
    p1.keep   	<- p2.keep <- abs(std.er[(2+sum(num.betas[1:5])):sum(num.betas[1:6])])<large
    if(het==T)	p2.keep   	<- abs(std.er[(2+sum(num.betas[1:6])):(sum(num.betas[1:6])+num.betas[6])])<large
    p1.cov		<- drop.bad.coef(p1.keep & p2.keep, p1.cov)
  }

}	# end the if statement



### drop covariates with the wrong sign (if no neg SE)
# before, I only dropped signif coefficients, but now drop all
if (mean(c(psi.keep, th0.keep, th1.keep, gam.keep, eps.keep, p1.keep, p2.keep, 1))==1) {

  if(num.betas[1]>1) {
    quads <- grep("psi.sq",betas)
    #psi.keep[seq(2, num.betas[1], by=2)] 	<- coefs[seq(3,num.betas[1], by=2)]/std.er[seq(3,num.betas[1], by=2)]<1.94
    psi.keep[quads-1] 	<- coefs[quads]<0	# '-1' is for intercept
    psi.cov	<- drop.bad.coef(psi.keep, psi.cov)
  }
  if(num.betas[4]>1) {
    quads <- grep("gam1.sq",betas)
    #gam.keep[seq(2, num.betas[4], by=2)] 	<- coefs[seq(sum(num.betas[1:3])+3, sum(num.betas[1:4]), by=2)]/std.er[seq(sum(num.betas[1:3])+3, sum(num.betas[1:4]), by=2)]<1.94
    gam.keep[quads-1-sum(num.betas[1:3])] 	<- coefs[quads]<0  # quads indicates position in all betas, but gam.keep is just for gam
    gam.cov	<- drop.bad.coef(gam.keep, gam.cov)
  }
  if(num.betas[5]>1) {
    quads <- grep("eps1.sq",betas)
    #eps.keep[seq(2, num.betas[5], by=2)] 	<- coefs[seq(sum(num.betas[1:4])+3, sum(num.betas[1:5]), by=2)]/std.er[seq(sum(num.betas[1:4])+3, sum(num.betas[1:5]), by=2)]>(-1.94)
    eps.keep[quads-1-sum(num.betas[1:4])] 	<- coefs[quads]>0
    eps.cov	<- drop.bad.coef(eps.keep, eps.cov)
  }

}	# end the if statement


initvals.index <- c(T,psi.keep, T,th0.keep, T,th1.keep, T,gam.keep, T,eps.keep, T,p1.keep)
if(het==T)  initvals.index <- c(initvals.index, T, p2.keep, T)
initvals  <- coefs[initvals.index]


### if no problem, finish analysis
if (mean(c(psi.keep, th0.keep, th1.keep, gam.keep, eps.keep, p1.keep, p2.keep, 1))==1) finished <- 1

### also check if intercepts have neg var
wonky2 <- is.na(sum(std.er))

# this checks for significant digits
jj=grep('significant digits', model.out)
num.conv <- model.out[jj]
sig.dig <- as.numeric(unlist(strsplit(num.conv, "\\s+"))[2])
wonky3 <- sig.dig<3     # TRUE if sig.dig low

print(paste0('cycle ', cycle))
cycle <- cycle+1

# if wonky result, move it to discard folder
if (finished==0 | wonky2==1 | wonky3==1) {
  file.rename(paste0('pres_',modname,'.out'), paste0('discard\\pres_',modname,'.out'))
}

if (cycle>40) finished <- 1

}
