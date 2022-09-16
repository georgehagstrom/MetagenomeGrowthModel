// Statistical model for the distribution of replication slopes for Synechococcus in Larkin et al  




data {

// In the data block we describe the types of data that this model takes as input


int<lower=0> NObsHLII;   // Number of observations of HLII metagenomes

int<lower=0> NDaysHLII;  // Number of unique days (measured noon-noon) on which those observations take place. 

// Note: Although each observation takes place in a different location, we approximate observations taking place on the same noon-noon day as occuring at the same place, and we fit
// a single replication slope function to the observations in that bin

vector[NObsHLII] RepSlopeVecHLII; // Calculated replication slope for each metagenome
vector[NObsHLII] RepSlopeSDVecHLII; // Standard deviation of the replication slope calculation function
vector<lower=0.0,upper=24.0>[NObsHLII] AdjustedTimeVecHLII; // This is the vector of the hour of each observation, measured noon-noon and based on local solar noon
int<lower=0> DayVecHLII[NObsHLII]; // This is the day of each observation


}


parameters {

// Here we define some of the parameters of the statistical model


real<lower=0.0> c0; // baseline slope- this the typical replication slope observed outside the time window for cell division


real<lower=0.0> MeanRepMax; //Replication slope max- this is 
real<lower=4.0,upper=10.0> t1; // Time where replication slope is maximal
real<lower=1.0,upper=6.0> tw; // The width of the replication interval
real <lower=0.0,upper=1.0> sigma; // the measurement variance
vector[NDaysHLII] MaxRepDayHLII; // The maximum replication slope on a given day 

}


transformed parameters {
vector[NObsHLII] ResponseFunctionHLII; // The response function is the modeled replication slope as a function of time on each day. 

vector[NObsHLII] RepSlopeHLII_star; // This is the difference between the replication slope and the response function

vector[NDaysHLII] muDayHLII; //Replication Slope at each site

for (n in 1:NObsHLII){  
    
// The functional form, we model the replication slope as a Gaussian centered on the time with maximum replication slope
    ResponseFunctionHLII[n] = c0 + (MaxRepDayHLII[DayVecHLII[n]+1]-c0)*exp(-(AdjustedTimeVecHLII[n]-t1)^2/(2.0*tw^2));

}

RepSlopeHLII_star = (ResponseFunctionHLII-RepSlopeVecHLII)/1;



}


model {

// start with priors

t1 ~ normal(10,2.0) T[0.0,24.0] ; // Time of maximal replication
c0 ~ normal(0.05,0.1); // Baseline replication slope

for (n in 1:NDaysHLII){
}

// We adopt a hierarchical model for the mean replication slope, assuming the distribution of the mean replication
// slope is a truncated normal. Priors are weakly informative

MeanRepMax ~ normal(0.2,1.0) T[0,];

// Each day has a modeled maximum replication slope, normally distributed around the mean.

MaxRepDayHLII ~ normal(MeanRepMax,0.3    ) ;

// RepSlopeHLII_star is the difference between the measured rep slope and the response function, which incorporates the
// diel cycle. sigma is the variance of this measurement

for (n in 1:NObsHLII){
RepSlopeHLII_star[n] ~ normal(0.0,sigma); 
}



}





// Generating these for posterior predictive checks

generated quantities {


vector[NObsHLII] RepSlopePostHLII;
vector[NDaysHLII] MaxRepSlopeHLII;




for (n in 1:NObsHLII){

RepSlopePostHLII[n] = normal_rng(ResponseFunctionHLII[n],sigma);

}


for (n in 1:NDaysHLII){

MaxRepSlopeHLII[n] = MaxRepDayHLII[n];


}


}






