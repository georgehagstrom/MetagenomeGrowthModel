
data {

int<lower=0> NObsHLII;
int<lower=0> NDaysHLII;
//vector<lower=0.0,upper=24.0>[NObsHLII] timeVecHLII;
vector[NObsHLII] RepSlopeVecHLII;
vector[NObsHLII] RepSlopeSDVecHLII;
vector<lower=0.0,upper=24.0>[NObsHLII] AdjustedTimeVecHLII;
int<lower=0> DayVecHLII[NObsHLII];


}


parameters {

real<lower=0.0> c0; // baseline slope


real<lower=0.0> MeanRepMax; //Replication slope max
real<lower=4.0,upper=10.0> t1; // Replication Onset
real<lower=1.0,upper=6.0> tw;
//real<lower=0.0,upper=1.0> sigmaM;
real <lower=0.0,upper=1.0> sigma;
//real<lower=0.0> sigmaRep;
vector[NDaysHLII] MaxRepDayHLII;
//vector<lower=-2,upper=2>[NDaysHLI] DeltaTVecHLI;

//vector<lower=-3,upper=3>[NDaysHLII] DeltaTVecHLII;
//vector<lower=-2,upper=2>[NDaysHLI] DeltaT1VecHLI;

//vector<lower=-3,upper=3>[NDaysHLII] DeltaT1VecHLII;




//vector<lower=-2,upper=2>[NDays] DeltaTVec;
//vector<lower=-2,upper=2>[NDays] DeltaT1Vec;


}


transformed parameters {

//vector[NObsHLI] RepSlopeModelHLI;

//vector[NObsHLII] RepSlopeModelHLII;
//vector[NObsHLI] ResponseFunctionHLI;

vector[NObsHLII] ResponseFunctionHLII;

//real t1_star;
//real c0_star;
//real muMean_star;

vector[NObsHLII] RepSlopeHLII_star;

vector[NDaysHLII] muDayHLII; //Replication Slope at each site
//real sigmaSite;
//real muMean = 0.7;
//real t1 = 18.0;
//real dt = 6.0;
//sigmaSite = 1.0;

//real sigma = 0.05;
//vector[NDays] muDay;
//real<lower=0.1,upper=0.4> sigmaSite;
//real<lower=0.2,upper=.4> sigma;

//sigmaSite = 1.0/InvSigmaSite;
//sigma = 1.0/InvSigma;
//vector<lower=-2,upper=2>[NDaysHLI] DeltaTVecHLI;

//vector<lower=-2,upper=2>[NDaysHLII] DeltaTVecHLII;
//vector<lower=-2,upper=2>[NDaysHLI] DeltaT1VecHLI;

//vector<lower=-2,upper=2>[NDaysHLII] DeltaT1VecHLII;


//t1_star = (t1-20.0)/4.0;
//c0_star = (c0-1.0);
//muDayHLII = (muDiffHLII_star*sigmaSiteHLII+muMean);



for (n in 1:NObsHLII){  
    

    ResponseFunctionHLII[n] = c0 + (MaxRepDayHLII[DayVecHLII[n]+1]-c0)*exp(-(AdjustedTimeVecHLII[n]-t1)^2/(2.0*tw^2));

 //   RepSlopeModelHLII[n] = log(exp(ResponseFunctionHLII[n]*k)+exp(k*c0))/k ;
}

RepSlopeHLII_star = (ResponseFunctionHLII-RepSlopeVecHLII)/1;



}


model {

// start with priors

//nu ~ normal(1.0,1.0) T[1,];
//sigma ~ exponential(1.0) T[0.1,];
//sigma ~ normal(0.1,0.1) T[0.1,3.0];
t1 ~ normal(10,2.0) T[0.0,24.0] ;
c0 ~ normal(0.05,0.1);
//dt ~ normal(6.0,1.0) T[3.0,8.0];

for (n in 1:NDaysHLII){
//DeltaTVecHLII[n] ~ normal(0.0,1.0) T[-3,3];
//DeltaT1VecHLII[n] ~ normal(0.0,1.0) T[-3,3];
}


//muDiff_star ~ normal(0.0,1.0);
//dt ~ normal(7.0,0.5) T[1,12];
//muDiffHLI_star ~ normal(0.0,1.0);

MeanRepMax ~ normal(0.2,1.0) T[0,];

MaxRepDayHLII ~ normal(MeanRepMax,0.3    ) ;

//muDiffHLII_star ~ normal(0.0,1.0);
//for (n in 1:NDays){

#muDay[n] ~ normal(muMean,0.3) ;

//muDay[n] ~ gamma(muMean^2/(sigmaSite)^2,muMean/(sigmaSite^2));

//} 

for (n in 1:NObsHLII){
RepSlopeHLII_star[n] ~ normal(0.0,sigma); //  student_t(5.0,0.0,1.0);
//RepSlope_star[n] ~ normal(0.0,1.0);
}



}







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






