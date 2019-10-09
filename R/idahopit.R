#program pit_stray
# C. Busack, NMFS
# May 11, 2017
# This program was developed to provide perspective on PIT detections of Idaho steelhead
# as an indicator of stray rates.  A preliminary summary of PIT detections over several years 
# yielded few detections.  One conclusion from a pattern of no or few detections over a 
# period of several years is that straying must be low, but given low tagging rates, the 
# the question is do few or no recoveries really mean few or no strays?
#
# This is a very simple simulation of what to expect in terms of PIT tag detections over a specified  
# number of years at one site, given that a specified number of returning adults are tagged at
# a specified rate and assuming 100% detection efficiency.  For example, if fish from a given group are 2% tagged, and 10 fish from 
# that group pass the detector each year, over a five year period, what is the probability of seeing 
# only no more than one detection during that period. Turns out that this probability is 
# about 74%

#------------------------------------------------------------------------------------
# Input
#------------------------------------------------------------------------------------
years <- 10
prob  <- 0.05  # tag rate of returning fish  
fish  <- 50    # number of fish from tagged group passing detector each year
#------------------------------------------------------------------------------------

trials <- 100000   # replicates; can change but no real point in doing so
w <- matrix(0, trials, years+1)

for (i in 1:years) {
  w[,i]<-rbinom(trials,fish,prob)  
}

for (i in 1:trials) {
  w[i,years+1]<-sum(w[i,1:years])
}

hist(w[,years + 1],xlab='Total Detections',main='PIT detections over years')

x <- table(w[,years+1])/trials
