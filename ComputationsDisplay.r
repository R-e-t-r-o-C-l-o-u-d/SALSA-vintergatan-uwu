if(!require(pracma)){
    install.packages("pracma")
}

if(!require(FITSio)){
    install.packages("FITSio")
}

if(!require(DescTools)){
    install.packages("DescTools")
}

library(pracma, help, pos = 2, lib.loc = NULL)
library(FITSio, help, pos = 2, lib.loc = NULL)
library(DescTools, help, pos = 2, lib.loc = NULL)

xplots <- c(0)
yplots <- c(0)

isEmpty <- function(x) {
    return(length(x)==0)
}


lookup <- function(hdr, keyword) {
    i <- which(hdr$key == keyword)
    if (length(i) == 0) {
        value <- NA
    } else {
        value <- hdr$value[i[1]]
    }
    value
}

filenames <- (Sys.glob("fits/*.fits"))
pb = txtProgressBar(min = 0, max = length(filenames), initial = 0, style=3, char="=", width=50) #progressbar
stepi=0

for(file in filenames){
stepi= stepi + 1
setTxtProgressBar(pb,stepi)

#reads fits file
f <- readFITS(file = file, hdu = 1, maxLines = 5000,
fixHdr = c('none', 'remove', 'substitute'), phdu = 1)


#defines hdr so it can be used for lookups
even <- seq(2, length(f$hdr), by=2)
odd <- even-1
hdr <- data.frame(key=f$hdr[odd], value=f$hdr[even], stringsAsFactors=FALSE) #defines dataframe for easy lookup of 

#defines values to calculate velocity for later
MHz <- 1.0e6              # define MHz in Hz
c.kms <- 2.99792458e5     # the speed of light in km/s
f0 <- as.double(lookup(hdr, "CRVAL1"))/MHz
df <- as.double(lookup(hdr, "CDELT1"))
vs <- as.double(lookup(hdr, "VLSR"))

if (is.na(vs)) vs <- as.double(lookup(hdr, "VELO-LSR"))
values <- c(f$imDat)


#look up for relevant values
dv <- -df/f0*c.kms/MHz
v <- dv*(seq(f$axDat$len[1])-f$axDat$crpix[1])-vs
vels <- c(0)
peaks <- (findpeaks(as.double(f$imDat), nups=2, ndown=4, minpeakheight=20,minpeakdistance=10,))
for(peak in peaks) {
    pe <- v[(which(values == peak))]
    if(!isEmpty(pe) && peak>0){
        vels <- append(vels,pe)
    }

    
}

l <- (as.double(lookup(hdr, "CRVAL2")))
conv<-(pi/180)
l <- l*conv
vels <- vels[-1]
R0 <-as.double(8.5)
V0 <-220
#defines constant values


for(vel in vels){

#geometry math, check paper if you want to understand
R <- (R0*V0*sin(l))/(V0*sin(l)+vel)

r <- sqrt(R**2-R0**2*sin(l)**2)+R0*cos(l)

#removes annoying stuff that makes program crash
if(is.nan(r)){
next()
}
xplots <- append(xplots, r*cos(l-pi/2))
yplots <- append(yplots, r*sin(l-pi/2)+8.5)
}

}
close(pb)

xplots <- xplots[-1]
yplots <- yplots[-1]
xplots <- append(xplots, 0)
yplots <- append(yplots, 8.5)
plot(yplots, xplots, xlab="Distance [kpc]",ylab="Distance [kpc]", xlim=c(-10,10), ylim=c(-0,20), lwd=1, pch=20, cex=0.455)
print("DONE")