#include <oxstd.h> 
#include <oxprob.h>

//Total Info;
decl totalN=2000;
decl NGRO=2;

//Group-level Info;
decl grp1N=1000, grp1J=10;
decl grp2N=1000, grp2J=10;
decl grp1data="Data1.dat";
decl grp2data="Data2.dat";
decl baseline=1;	// 1=Constrained-baseline
					// 2=Free-baseline
decl totalJ=10;		// Change total number of items based on the baseline model
					// If Constrained-baseline, totalJ = grp1J
					// If Free-baseline, totalJ = grp1J+grp2J-anchor
					
decl anchor={1,2};	// If free-baseline is used, indicate anchor items
decl test=3;		// If free-baseline is used, indicate test item

decl time0, time1;
decl X;
decl theta0,alpha0,beta0,mn0,sd0;
decl Alphaa=1.5, Alphab=1.5, Alphamin=0.25, Alphamax=4;	 // 4-parameter beta prior distribution parameters for alpha
decl Betaa=2, 	 Betab=2, 	 Betamin=-5, 	Betamax=5;	 // 4-parameter beta prior distribution parameters for beta
decl cand_th=.90,cand_a=.20,cand_b=.20;					 // proposal distribution standard deviations
decl cand_mn=.02;  //The focal group subpopulation mean candidate sampling distribution standard deviation
decl cand_sd=.02;  //The focal group subpopulation SD candidate sampling distribution standard deviation
decl acc,accind;
decl acpttheta,acpta,acptb,acptmn,acptsd;
decl Grp; //N-size vector of group membership;
decl maxsecs=999999;	  //Maximum # of seconds for estimation to run. Note: 86400 is 24 hours;
decl maxmins=720; //Maximum # of minutes for estimation to run;
decl Counter=500;		//After how many updates do you want to receive notification?
decl it=10000,burn=5000;
decl chains=3;
decl output=0;

decl EstMean=1;	//Allow the population means to vary (Group 1 is fixed at N(0,1);
decl EstSD=1;	//Allow the population standard deviations to vary (Group 1 is fixed at N(0,1));

decl outputresults=0; //1=output item parameter results to file, 0=Do not output;
	decl resultfile="ParEsts.xls";//File name for parameter results
decl outputtheta=0; //1=output theta estimates to file, 0=Do not output;
	decl thetafile="ThetaEsts.xls";  //File name for estimated thetas

decl outputsamples=1; //1=output MC samples, 0=Do not output;
decl outputINsamplesG1="MCSamples_infit_G1.txt";
decl outputOUTsamplesG1="MCSamples_outfit_G1.txt";
decl outputOSDsamplesG1="MCSamples_osd_G1.txt";
decl outputQ1samplesG1="MCSamples_q1_G1.txt";
decl outputINsamplesG2="MCSamples_infit_G2.txt";
decl outputOUTsamplesG2="MCSamples_outfit_G2.txt";
decl outputOSDsamplesG2="MCSamples_osd_G2.txt";
decl outputQ1samplesG2="MCSamples_q1_G2.txt";

//======= Files to read priors from =======
decl AlphaPriorFile="PriorsAlpha.txt";
decl BetaPriorFile="PriorsBeta.txt";

//==========  Begin program functions =============;

TwoPL(const a,b,th){
	decl sN = rows(th);
	return((1)./(1+exp(-a'.*(th-ones(sN,1)*b'))));
}		 

like(const a,b,th){
	decl P=TwoPL(a,b,th);
	decl tmpLike=X.*log(P)+(1-X).*log((1-P));
	decl Like = isdotnan(tmpLike) .? 0.0 .: tmpLike; // replace NaN with 0.0
	return(Like);
}


initiate(const c){

	//THETA INTIALS from RandN(0,1) for thetas;
	theta0=rann(totalN,1);  

	//PARAMETER INITIALS
	alpha0=ones(totalJ,1);
	beta0=zeros(totalJ,1);

	decl inits=zeros(totalJ,2);
	for(decl p=0;p<totalJ;p++){
		inits[p][0]=alpha0[p];
		inits[p][1]=beta0[p];
	}

	mn0=zeros(1,NGRO);
	sd0=ones(1,NGRO);
	
	if(c==0){println("Initial Parameter Values:"); print(inits);}
}

drawtheta(){
	decl d0,d1;
	decl theta1=cand_th*rann(totalN,1)+theta0;
	decl thmu=zeros(totalN,1);
	decl thSD=ones(totalN,1);
	for(decl p=0;p<totalN;p++){
		thmu[p]=mn0[Grp[p]];
		thSD[p]=sd0[Grp[p]];
	}

	d0=sumr(like(alpha0,beta0,theta0))+log(densn((theta0-thmu)./thSD));	
	d1=sumr(like(alpha0,beta0,theta1))+log(densn((theta1-thmu)./thSD));	
	acc=d1-d0;
	accind=vecindex(exp(acc).>ranu(totalN,1));
	theta0[accind]=theta1[accind];		   //Only updates indices in the vector where exp(d1-d0) > a random uniform;
	acpttheta[accind]=acpttheta[accind]+1; //If value accepted, adds one to the tracking acceptance rate;
}

drawbeta(){
	decl d0,d1;
	decl alpha1=cand_a*rann(totalJ,1)+alpha0;  	 
	d0=sumc(like(alpha0,beta0,theta0))+log(densbeta((alpha0'-Alphamin')./(Alphamax'-Alphamin'),Alphaa',Alphab'));
	d1=sumc(like(alpha1,beta0,theta0))+log(densbeta((alpha1'-Alphamin')./(Alphamax'-Alphamin'),Alphaa',Alphab'));
	acc=d1-d0;
	accind=vecindex(exp(acc).>ranu(1,totalJ));
	alpha0[accind]=alpha1[accind];
	acpta[accind]=acpta[accind]+1;		 //If value accepted, adds one to the tracking acceptance rate;

	decl d2,d3;
	decl beta1=cand_b*rann(totalJ,1)+beta0;	
	d2=sumc(like(alpha0,beta0,theta0))+log(densbeta((beta0'-Betamin')./(Betamax'-Betamin'),Betaa',Betab'));
	d3=sumc(like(alpha0,beta1,theta0))+log(densbeta((beta1'-Betamin')./(Betamax'-Betamin'),Betaa',Betab'));
	acc=d3-d2;		
	accind=vecindex(exp(acc).>ranu(1,totalJ));
	beta0[accind]=beta1[accind];
	acptb[accind]=acptb[accind]+1;
}



//Procedure to draw and population parameters (mean, SD);
drawlamda(){
	decl d0,d1;

	//Populate two vectors of that theta's group mean & SD;
	decl thmu0=zeros(totalN,1);
	decl thSD0=ones(totalN,1);
	for(decl p=0;p<totalN;p++){
		thmu0[p]=mn0[Grp[p]];
		thSD0[p]=sd0[Grp[p]];
	}

	//Draw mean candidates;
	decl mn1=cand_mn*rann(1,NGRO)+mn0;
	decl sd1=cand_sd*rann(1,NGRO)+sd0;
		//Group 1 set to N(0,1);
		mn1[0]=0; sd1[0]=1;
	decl thmu1=zeros(totalN,1);
	decl thSD1=ones(totalN,1);
	for(decl p=0;p<totalN;p++){
		thmu1[p]=mn1[Grp[p]];
		thSD1[p]=sd1[Grp[p]];
	}

	//------- First, test means ------
	//If user selected, test mean candidates;
	if(EstMean==1){
		//like is P(thetas | Mn,SD), uniform prior results in the prior dropping out;
		d0=sumc(log(densn((theta0-thmu0)./thSD0)./thSD0));
		d1=sumc(log(densn((theta0-thmu1)./thSD0)./thSD0));

		//Update means based on full like;
		acc=d1-d0;
		
		if(exp(acc).>ranu(1,totalJ)){	//Updates new means where exp(model1-model0) > a random uniform;
			mn0=mn1;
			acptmn=acptmn+1;
		}
	} //End EstMean If;

	//------- Second, test SDs ------;
	//If user selected, test SD candidates;
	if(EstSD==1){
		//like is P(thetas | Mn,SD), uniform prior results in the prior dropping out;
		d0=sumc(log(densn((theta0-thmu0)./thSD0)./thSD0));
		d1=sumc(log(densn((theta0-thmu0)./thSD1)./thSD1));

		//Update means based on full like;
		acc=d1-d0;

		
		if(exp(acc).>ranu(1,totalJ)){	//Updates new means where exp(model1-model0) > a random uniform;
			sd0=sd1;
			acptsd=acptsd+1;
		}
	} //End EstSD If;
}

//========================================Begin Main Program=========================================;
main()                
{
time0 = timer();
ranseed(today());

//Resize matrices;
Grp=zeros(totalN,1);
X=constant(.NaN,totalN,totalJ);

//====Read in group data====;
decl grp1x, grp2x; //Temporary matrices to store group-level data;
//---Group 1---
decl g1file=fopen(grp1data);
fscan(g1file,"%#M",grp1N,grp1J,&grp1x);
fclose(g1file);
//---Group 2---
decl g2file=fopen(grp2data);
fscan(g2file,"%#M",grp2N,grp2J,&grp2x);
fclose(g2file);

//Place group membership into Grp vector;
Grp[0:grp1N-1]=0;
Grp[grp1N:grp1N+grp2N-1]=1;

//---Input items related to each group---
//decl grp1its,grp2its;
//decl gitfile=fopen(grpitems);
//fscan(gitfile,"%#M",1,grp1J,&grp1its);
//fscan(gitfile,"%#M",1,grp2J,&grp2its);
//fclose(gitfile);

decl grp1its=zeros(1,grp1J);
decl grp2its=zeros(1,grp1J);
for(decl i=0; i<grp1J; i++){
	grp1its[i]=i+1;
}
decl grpcounter=1;
if(baseline==1){
	grp2its=grp1its;
}
if(baseline==2){
	for(decl i=0; i<grp2J; i++){
		if(i+1==anchor[0] || i+1==anchor[1] || i+1==test){
			grp2its[i]=i+1;
			grpcounter=grpcounter;
		}
		else{
			grp2its[i]=grp1J+grpcounter;
			grpcounter=grpcounter+1;
		}
	}
}

//PLace response data for the groups into a single X response data matrix based on which items they had;
X[0:grp1N-1][grp1its-1]=grp1x[][0:grp1J-1];
X[grp1N:grp1N+grp2N-1][grp2its-1]=grp2x[][0:grp2J-1];

//Output first and last lines of data;
println("Responses for first line of data");
print(X[0][]);
println("Responses for last line of data");
print(X[totalN-1][]);

//Output Alpha priors to screeen;
Alphaa=Alphaa*ones(totalJ,1);
Alphab=Alphab*ones(totalJ,1);
Alphamin=Alphamin*ones(totalJ,1);
Alphamax=Alphamax*ones(totalJ,1);
println("Alpha priors");
print( "%c", {"Item", "Shape a", "Shape b", "Min val", "Max val"}, "%cf",{"%8.4g", "%8.4g", "%8.4g", "%8.4g", "%8.4g"}, range(1, totalJ, 1)'~Alphaa~Alphab~Alphamin~Alphamax);

Betaa=Betaa*ones(totalJ,1);
Betab=Betab*ones(totalJ,1);
Betamin=Betamin*ones(totalJ,1);
Betamax=Betamax*ones(totalJ,1);
//Output Delta priors to screeen;
println("Beta priors");
print( "%c", {"Item", "Shape a", "Shape b", "Min val", "Max val"}, "%cf",{"%8.4g", "%8.4g", "%8.4g", "%8.4g", "%8.4g"}, range(1, totalJ, 1)'~Betaa~Betab~Betamin~Betamax);


//Matrices to store the state of each chain (elements set to zero) during sampling;
decl thstate=zeros(totalN,chains);
decl astate=zeros(totalJ,chains);
decl bstate=zeros(totalJ,chains);
decl mnstate=zeros(chains,NGRO);
decl sdstate=zeros(chains,NGRO);

//Matrices & Variables to store a running sum for each chain (elements set to zero);
decl sthchains=zeros(totalN,chains);
decl sachains=zeros(totalJ,chains);
decl sbchains=zeros(totalJ,chains);
decl smnchains=zeros(chains,NGRO);
decl ssdchains=zeros(chains,NGRO);
decl sthchains2=zeros(totalN,chains);
decl sachains2=zeros(totalJ,chains);
decl sbchains2=zeros(totalJ,chains);
decl smnchains2=zeros(chains,NGRO);
decl ssdchains2=zeros(chains,NGRO);
decl sthsum,sasum,sbsum;
decl sthsum2,sasum2,sbsum2;
decl smnsum,ssdsum;
decl smnsum2,ssdsum2;
acpta=zeros(totalJ,1);
acptb=zeros(totalJ,1);
acpttheta=zeros(totalN,1);
acptmn=zeros(NGRO,1);
acptsd=zeros(NGRO,1);

//Matrices to store complete run of chains;
decl alphastates=zeros(totalJ*chains,it-burn);
decl betastates=zeros(totalJ*chains,it-burn);
decl thetastates=zeros(totalN*chains,it-burn);

//Initiate state 0 for each variable & chain;
for(decl c=0;c<chains;c++){
	//Obtain initial values for each variable from the function;
	initiate(c);
	thstate[][c]=theta0;
	astate[][c]=alpha0;
	bstate[][c]=beta0;
	mnstate[c][]=mn0;
	sdstate[c][]=sd0;
}

//Variables to track iterations and flag whether to continue;
decl Cont=1;
decl i=0;

//Begin iterations for the chains;
while(Cont==1) {
	//To provide feedback during update: Print every Counter update;
	if(imod(i+1,Counter)==0) println("Beginning iteration ", i+1, ". Lapsed time so far: ", timespan(time0));
	for(decl c=0;c<chains;c++){
		//Set the temporary vectors to the current state of the chain;
		mn0=mnstate[c][];
		sd0=sdstate[c][];
		theta0=thstate[][c];
		alpha0=astate[][c];
		beta0=bstate[][c];
		//Now draw & update;
		drawtheta();
		drawbeta();		  	
		drawlamda();
		if(i>burn-1){  //Note: v4 removed call to update function, updates here;
			sthchains[][c]=sthchains[][c]+theta0;
			sthchains2[][c]=sthchains2[][c]+theta0.^2;
			sachains[][c]=sachains[][c]+alpha0;
			sachains2[][c]=sachains2[][c]+alpha0.^2;
			sbchains[][c]=sbchains[][c]+beta0;
			sbchains2[][c]=sbchains2[][c]+beta0.^2;
			smnchains[c][]=smnchains[c][]+mn0;
			smnchains2[c][]=smnchains2[c][]+mn0.^2;
			ssdchains[c][]=ssdchains[c][]+sd0;
			ssdchains2[c][]=ssdchains2[c][]+sd0.^2;

			//Save state of chain for all iterations in full matrix for later use;
			alphastates[c*totalJ:c*totalJ+(totalJ-1)][i-burn]=alpha0;
			betastates[c*totalJ:c*totalJ+(totalJ-1)][i-burn]=beta0;
			thetastates[c*totalN:c*totalN+(totalN-1)][i-burn]=theta0;
		}
		
		//Now, save the current state for the chain for next cycle;
		thstate[][c]=theta0;
		astate[][c]=alpha0;
		bstate[][c]=beta0;
		mnstate[c][]=mn0;
		sdstate[c][]=sd0;
		
	} //End chain, update the next chain;

	//Functions to determine if updating should continue;
	i=i+1;
	if(i==it) Cont=0;
	time1=timer();
	//if(time1-time0>maxmins*100*60) Cont=0;
} //End iteration loop;

//////////////////////////////////////////////////////////////////////////////////////////
/*
//Output results to screen;
println("Number of Iterations Performed: ", i);
for(decl c=0;c<chains;c++){
	println("Chain ", c+1, " Alpha Parameter Means:");
	print( "%c", {"Alpha Parameter ", "Mean"}, "%cf",{"%8.4g", "%#13.5g"}, range(1, totalJ, 1)'~(sachains[][c]/(i-burn)));
	println("Chain ", c+1, " Beta Parameter Means:");
	print( "%c", {"Beta Parameter ", "Mean"}, "%cf",{"%8.4g", "%#13.5g"}, range(1, totalJ, 1)'~(sbchains[][c]/(i-burn)));
		  
	println("Chain ", c+1, " Parameter Alpha SDs:");
	print( "%c", {"Alpha Parameter", "SD"}, "%cf",{"%8.4g", "%#13.5g"}, range(1, totalJ, 1)'~(sqrt(((sachains2[][c])-((sachains[][c]).^2)/(i-burn))/(i-burn-1))));
	println("Chain ", c+1, " Parameter Delta SDs:");
	print( "%c", {"Beta Parameter", "SD"}, "%cf",{"%8.4g", "%#13.5g"}, range(1, totalJ, 1)'~(sqrt(((sbchains2[][c])-((sbchains[][c]).^2)/(i-burn))/(i-burn-1))));
	
	println("Chain ", c+1, " Group Mean & Group SD");
	print( "%c", {"Group", "Mean", "SD"}, "%cf",{"%8.4g", "%#13.5g", "%#13.5g"},  range(0, NGRO-1, 1)'~((smnchains[c][])/(it-burn))'~((ssdchains[c][])/(it-burn))');
	println("Chain ", c+1, " SD for the group parameters (Mean & SD)");
	print( "%c", {"Group", "SD (Mean)", "SD (SD)"}, "%cf",{"%8.4g", "%#13.5g", "%#13.5g"},  range(0, NGRO-1, 1)'~(sqrt(((smnchains2[c][])-((smnchains[c][]).^2)/(i-burn))/(i-burn-1)))'~(sqrt(((ssdchains2[c][])-((ssdchains[c][]).^2)/(i-burn))/(i-burn-1)))');
}
*/

//Sum across the chains, for calculating mean & SD across chains
sthsum=sumr(sthchains);
sasum=sumr(sachains);
sbsum=sumr(sbchains);
sthsum2=sumr(sthchains2);
sasum2=sumr(sachains2);
sbsum2=sumr(sbchains2);
smnsum=sumc(smnchains);
ssdsum=sumc(ssdchains);
smnsum2=sumc(smnchains2);
ssdsum2=sumc(ssdchains2);

//--- Output results to screen ---;
println("------ Alpha Means, side-by-side for each chain ------");
println(sachains/(i-burn));
println("------ Beta Means, side-by-side for each chain ------");
println(sbchains/(i-burn));
println("------ Mean across chains (Alpha Beta) ------");
decl means=(sasum~sbsum)/((i-burn)*chains);
println(means);
println("------ Standard Deviation across chains (Alpha Beta) ------");
decl SDs=sqrt(((sasum2~sbsum2)-((sasum~sbsum).^2)/((i-burn)*chains))/((i-burn-1)*chains));
println(SDs);
println("------ Group Means across chains  ------");
print( "%c", {"Group", "Mean", "SD"}, "%cf",{"%8.4g", "%#13.5g", "%#13.5g"},  range(0, NGRO-1, 1)'~((smnsum)/((it-burn)*chains))'~((ssdsum)/((it-burn)*chains))');


////////////////////Calculate R statistics//////////////////

decl Ba=zeros(totalJ,chains);
decl Bb=zeros(totalJ,chains);

for(decl c=0;c<chains;c++){
	Ba[][c]= ((sachains[][c]/(i-burn))-(sasum/((i-burn)*chains))).^2;
	Bb[][c]= ((sbchains[][c]/(i-burn))-(sbsum/((i-burn)*chains))).^2;
}
	
decl Balpha=sumr(Ba)/(chains-1);
decl Bbeta=sumr(Bb)/(chains-1);

decl Wa=zeros(totalJ*chains,it-burn);
decl Wb=zeros(totalJ*chains,it-burn);

for(decl c=0;c<chains;c++){
	Wa[c*totalJ:c*totalJ+(totalJ-1)][]=(alphastates[c*totalJ:c*totalJ+(totalJ-1)][]-sachains[][c]/(i-burn)).^2;
	Wb[c*totalJ:c*totalJ+(totalJ-1)][]=(betastates[c*totalJ:c*totalJ+(totalJ-1)][]-sbchains[][c]/(i-burn)).^2;
}

decl Waa=sumr(Wa);
decl Wbb=sumr(Wb);

decl Waaa=zeros(totalJ,chains);
decl Wbbb=zeros(totalJ,chains);

for(decl c=0;c<chains;c++){
	Waaa[][c]=Waa[c*totalJ:c*totalJ+(totalJ-1)];
	Wbbb[][c]=Wbb[c*totalJ:c*totalJ+(totalJ-1)];
}	
	
decl Walpha=sumr(Waaa)/((i-burn)*chains);
decl Wbeta=sumr(Wbbb)/((i-burn)*chains);

decl Salpha=(i-burn-1)/(i-burn).*Walpha+Balpha;
decl Sbeta=(i-burn-1)/(i-burn).*Wbeta+Bbeta;

decl Ralpha=((chains+1)/chains).*Salpha./Walpha-((i-burn-1)/(chains*(i-burn)));
decl Rbeta=((chains+1)/chains).*Sbeta./Wbeta-((i-burn-1)/(chains*(i-burn)));

println("------ Convergence Criterion: R statistics < 1.2 ------");
println("Alpha Parameters:");
println(Ralpha);
println("Beta Parameters:");
println(Rbeta);

println("------ Acceptance Rates  ------");
println("Alpha Parameters:");
println((acpta)/(i*chains));
println("Beta Parameters:");
println((acptb)/(i*chains));
println("Group Means:");
println((acptmn)/(i*chains));
println("Group SDs:");
println((acptsd)/(i*chains));
println("Average theta:");
println(sumc(acpttheta/(i*chains))/totalN);

// Compute Likelihood Ratio //

decl alpha=means[][0];
decl beta=means[][1];
decl theta=sthsum/((i-burn)*chains);

decl LL = sumr(sumc(like(alpha,beta,theta)));
decl npar = 2*totalJ;
decl AIC = -2*LL+npar;
decl BIC = -2*LL+npar*log(totalN);
decl Dbar=0;

decl sumchainsalpha=zeros(totalJ,it-burn);
decl sumchainsbeta=zeros(totalJ,it-burn);
decl sumchainstheta=zeros(totalN,it-burn);

for(decl c=0;c<chains;c++){
	sumchainsalpha=sumchainsalpha+alphastates[(c*totalJ):(c*totalJ+(totalJ-1))][];
	sumchainsbeta=sumchainsbeta+betastates[(c*totalJ):(c*totalJ+(totalJ-1))][];
   	sumchainstheta=sumchainstheta+thetastates[(c*totalN):(c*totalN+(totalN-1))][]; 
}
	
decl avgchalpha=sumchainsalpha/chains;
decl avgchbeta=sumchainsbeta/chains;
decl avgchth=sumchainstheta/chains;

for (decl i=0;i<it-burn;i++){
	decl Dbartmp=sumr(sumc(like(avgchalpha[][i],avgchbeta[][i],avgchth[][i])));
	Dbar=Dbar+Dbartmp;
}
Dbar=-2*Dbar/(it-burn);
decl pD=Dbar-(-2*LL);
decl DIC=Dbar+pD;	

print("\n","-------        Model Fit Statistics        ------");
print("%r",{" log-Likelihood    "," Deviance    "," AIC    "," BIC    "," DIC    "},"%10.4f", LL|-2*LL|AIC|BIC|DIC);
print("\n");


////////////// Posterior Predictive Model Check //////////////

decl OSDpvalG1=zeros(grp1J,1);
decl OSDpvalG2=zeros(grp2J,1);
decl OSDsampleG1=zeros(grp1J,chains*(it-burn));
decl OSDsampleG2=zeros(grp2J,chains*(it-burn));

decl OUTpvalG1=zeros(grp1J,1);
decl OUTpvalG2=zeros(grp2J,1);
decl OUTsampleG1=zeros(grp1J,chains*(it-burn));
decl OUTsampleG2=zeros(grp2J,chains*(it-burn));

decl INpvalG1=zeros(grp1J,1);
decl INpvalG2=zeros(grp2J,1);
decl INsampleG1=zeros(grp1J,chains*(it-burn));
decl INsampleG2=zeros(grp2J,chains*(it-burn));

decl Q1pvalG1=zeros(grp1J,1);
decl Q1pvalG2=zeros(grp2J,1);
decl Q1sampleG1=zeros(grp1J,chains*(it-burn));
decl Q1sampleG2=zeros(grp2J,chains*(it-burn));

decl ChisqpvalG1=zeros(grp1J,1);
decl ChisqpvalG2=zeros(grp2J,1);
decl ChisqsampleG1=zeros(grp1J,chains*(it-burn));
decl ChisqsampleG2=zeros(grp2J,chains*(it-burn));

decl grid=zeros(61,1);
for(decl i=0; i<61; i++){
	grid[i]=-3+.1*i;
}
decl pgrid=zeros(61,1);
for(decl j=0; j<61; j++){
	pgrid[j]=densn(grid[j]);
}
decl counter=1;
for(decl c=0; c<chains; c++){
	for(decl rep=0; rep<it-burn; rep++){

		// Replicated Data Generation
		decl repalpha = alphastates[c*totalJ:c*totalJ+(totalJ-1)][rep];
		decl repbeta = betastates[c*totalJ:c*totalJ+(totalJ-1)][rep];
		decl repth = thetastates[c*totalN:c*totalN+(totalN-1)][rep];
		decl repaG1 = repalpha[grp1its-1];
		decl repaG2 = repalpha[grp2its-1];
		decl repbG1 = repbeta[grp1its-1];
		decl repbG2 = repbeta[grp2its-1];
		decl repthG1 = repth[0:grp1N-1];
		decl repthG2 = repth[grp1N:totalN-1];

		////////////////// Group 1 //////////////////
		// Data Generation
		decl repXG1=zeros(grp1N,grp1J);
		decl repPG1=TwoPL(repaG1,repbG1,repthG1);
		decl randG1=ranu(grp1N,grp1J);
		for(decl p=0;p<grp1N;p++){
			for(decl j=0;j<grp1J;j++){
				if(repPG1[p][j] > randG1[p][j]){ repXG1[p][j]=1; } 
				else{ repXG1[p][j]=0; }
			}
		}

		// Observed Discrepency Statistics
		
			// OSD
		decl G1obsNC=sumc(grp1x);
		decl G1obsENC=sumc(repPG1);
		decl G1obsOSD=((G1obsNC-G1obsENC).^2)./G1obsENC;

			// Outfit
		decl G1obsoutfitE = repPG1;
		decl G1obsoutfitp0 = 1-G1obsoutfitE;
		decl G1obsoutfitp1 = G1obsoutfitE;
		decl G1obsoutfitVar = (((0 - G1obsoutfitE).^2).*G1obsoutfitp0) + (((1 - G1obsoutfitE).^2).*G1obsoutfitp1);
		decl G1obsoutfitZ = (grp1x - G1obsoutfitE)./sqrt(G1obsoutfitVar);
		decl G1obsoutfitsq = G1obsoutfitZ.^2;
		decl G1obsoutfit = sumc(G1obsoutfitsq)/grp1N;

			// Infit
		decl G1obsinfit = sumc(G1obsoutfitVar.*G1obsoutfitsq)./sumc(G1obsoutfitVar);

			// Q1
		decl G1obstempmatrix=grp1x~repthG1;
		decl G1obsordermatrix=sortbyc(G1obstempmatrix,grp1J);
		decl G1obsorderX=dropc(G1obsordermatrix,grp1J);
		decl G1obsorderth=G1obsordermatrix[][grp1J];
		decl G1obsorderP=TwoPL(repaG1,repbG1,G1obsorderth);
		decl G1obsNN=grp1N/10;
		decl G1obsYensQ=zeros(1,grp1J);
		for(decl i=0;i<10;i++){
			decl G1obsbigO=(1/G1obsNN)*sumc(G1obsorderX[(G1obsNN*i):(G1obsNN*i+(G1obsNN-1))][]);
			decl G1obsbigE=(1/G1obsNN)*sumc(G1obsorderP[(G1obsNN*i):(G1obsNN*i+(G1obsNN-1))][]);
			G1obsYensQ=G1obsYensQ+G1obsNN*(((G1obsbigO-G1obsbigE).^2)./((G1obsbigE).*(1-G1obsbigE)));
		}

			// Chisquare
		decl G1obsP1out=TwoPL(repaG1,repbG1,grid); 
		decl G1obsP0out=1-G1obsP1out;
		decl G1obsPout=G1obsP0out~G1obsP1out;
		decl G1obsPP0out=grp1N*sumc(G1obsP0out.*(pgrid*ones(1,grp1J)))/10;
		decl G1obsPP1out=grp1N*sumc(G1obsP1out.*(pgrid*ones(1,grp1J)))/10;
		decl G1obsExp0=G1obsPP0out';
		decl G1obsExp1=G1obsPP1out';
		decl G1obsObs0=sumc((1-grp1x))';
		decl G1obsObs1=sumc(grp1x)';
		decl G1obsChisq0=((G1obsObs0-G1obsExp0).^2)./G1obsExp0;
		decl G1obsChisq1=((G1obsObs1-G1obsExp1).^2)./G1obsExp1;
		decl G1obsChisq=G1obsChisq0+G1obsChisq1;

		// Replicated Discrepency Statistics
		
			// OSD
		decl G1repNC=sumc(repXG1);
		decl G1repENC=sumc(repPG1);
		decl G1repOSD=((G1repNC-G1repENC).^2)./G1repENC;
		for(decl k=0;k<grp1J;k++){
			OSDsampleG1[k][counter-1]=G1repOSD[k];
		}
			// Outfit
		decl G1repoutfitE = repPG1;
		decl G1repoutfitp0 = 1-G1repoutfitE;
		decl G1repoutfitp1 = G1repoutfitE;
		decl G1repoutfitVar = (((0 - G1repoutfitE).^2).*G1repoutfitp0) + (((1 - G1repoutfitE).^2).*G1repoutfitp1);
		decl G1repoutfitZ = (repXG1 - G1repoutfitE)./sqrt(G1repoutfitVar);
		decl G1repoutfitsq = G1repoutfitZ.^2;
		decl G1repoutfit = sumc(G1repoutfitsq)/grp1N;
		for(decl k=0;k<grp1J;k++){
			OUTsampleG1[k][counter-1]=G1repoutfit[k];
		}
			// Infit
		decl G1repinfit = sumc(G1repoutfitVar.*G1repoutfitsq)./sumc(G1repoutfitVar);
		for(decl k=0;k<grp1J;k++){
			INsampleG1[k][counter-1]=G1repinfit[k];
		}
			// Yen's Q1
		decl G1reptempmatrix=repXG1~repthG1;
		decl G1repordermatrix=sortbyc(G1reptempmatrix,grp1J);
		decl G1reporderX=dropc(G1repordermatrix,grp1J);
		decl G1reporderth=G1repordermatrix[][grp1J];
		decl G1reporderP=TwoPL(repaG1,repbG1,G1reporderth);
		decl G1repNN=grp1N/10;
		decl G1repYensQ=zeros(1,grp1J);
		for(decl i=0;i<10;i++){
			decl G1repbigO=(1/G1repNN)*sumc(G1reporderX[(G1repNN*i):(G1repNN*i+(G1repNN-1))][]);
			decl G1repbigE=(1/G1repNN)*sumc(G1reporderP[(G1repNN*i):(G1repNN*i+(G1repNN-1))][]);
			G1repYensQ=G1repYensQ+G1repNN*(((G1repbigO-G1repbigE).^2)./((G1repbigE).*(1-G1repbigE)));
		}
		for(decl k=0;k<grp1J;k++){
			Q1sampleG1[k][counter-1]=G1repYensQ[k];
		}
		
			// Chisquare
		decl G1repP1out=TwoPL(repaG1,repbG1,grid);
		decl G1repP0out=1-G1repP1out;
		decl G1repPout=G1repP0out~G1repP1out;
		decl G1repPP0out=grp1N*sumc(G1repP0out.*(pgrid*ones(1,grp1J)))/10;
		decl G1repPP1out=grp1N*sumc(G1repP1out.*(pgrid*ones(1,grp1J)))/10;
		decl G1repExp0=G1repPP0out';
		decl G1repExp1=G1repPP1out';
		decl G1repObs0=sumc((1-repXG1))';
		decl G1repObs1=sumc(repXG1)';
		decl G1repChisq0=((G1repObs0-G1repExp0).^2)./G1repExp0;
		decl G1repChisq1=((G1repObs1-G1repExp1).^2)./G1repExp1;
		decl G1repChisq=G1repChisq0+G1repChisq1;
		for(decl k=0;k<grp1J;k++){
			ChisqsampleG1[k][counter-1]=G1repChisq[k];
		}

		// Compute the predictive posterior p-value (ppp)
		decl G1OSDind=vecindex(G1repOSD.>G1obsOSD);
		decl G1outfitind=vecindex(G1repoutfit.>G1obsoutfit);
		decl G1infitind=vecindex(G1repinfit.>G1obsinfit);
		decl G1YensQind=vecindex(G1repYensQ.>G1obsYensQ);
		decl G1Chisqind=vecindex(G1repChisq.>G1obsChisq);
		OSDpvalG1[G1OSDind]=OSDpvalG1[G1OSDind]+1;
		OUTpvalG1[G1outfitind]=OUTpvalG1[G1outfitind]+1;
		INpvalG1[G1infitind]=INpvalG1[G1infitind]+1;
		Q1pvalG1[G1YensQind]=Q1pvalG1[G1YensQind]+1;
		ChisqpvalG1[G1Chisqind]=ChisqpvalG1[G1Chisqind]+1;


		////////////////// Group 2 //////////////////
		// Data Generation
		decl repXG2=zeros(grp2N,grp2J);
		decl repPG2=TwoPL(repaG2,repbG2,repthG2);
		decl randG2=ranu(grp2N,grp2J);
		for(decl p=0;p<grp2N;p++){
			for(decl j=0;j<grp2J;j++){
				if(repPG2[p][j] > randG2[p][j]){ repXG2[p][j]=1; } 
				else{ repXG2[p][j]=0; }
			}
		}

		// Observed Discrepency Statistics

			// OSD
		decl G2obsNC=sumc(grp2x);
		decl G2obsENC=sumc(repPG2);
		decl G2obsOSD=((G2obsNC-G2obsENC).^2)./G2obsENC;

			// Outfit
		decl G2obsoutfitE = repPG2;
		decl G2obsoutfitp0 = 1-G2obsoutfitE;
		decl G2obsoutfitp1 = G2obsoutfitE;
		decl G2obsoutfitVar = (((0 - G2obsoutfitE).^2).*G2obsoutfitp0) + (((1 - G2obsoutfitE).^2).*G2obsoutfitp1);
		decl G2obsoutfitZ = (grp2x - G2obsoutfitE)./sqrt(G2obsoutfitVar);
		decl G2obsoutfitsq = G2obsoutfitZ.^2;
		decl G2obsoutfit = sumc(G2obsoutfitsq)/grp2N;

			// Infit
		decl G2obsinfit = sumc(G2obsoutfitVar.*G2obsoutfitsq)./sumc(G2obsoutfitVar);

			// Q1
		decl G2obstempmatrix=grp2x~repthG2;
		decl G2obsordermatrix=sortbyc(G2obstempmatrix,grp2J);
		decl G2obsorderX=dropc(G2obsordermatrix,grp2J);
		decl G2obsorderth=G2obsordermatrix[][grp2J];
		decl G2obsorderP=TwoPL(repaG2,repbG2,G2obsorderth);
		decl G2obsNN=grp2N/10;
		decl G2obsYensQ=zeros(1,grp2J);
		for(decl i=0;i<10;i++){
			decl G2obsbigO=(1/G2obsNN)*sumc(G2obsorderX[(G2obsNN*i):(G2obsNN*i+(G2obsNN-1))][]);
			decl G2obsbigE=(1/G2obsNN)*sumc(G2obsorderP[(G2obsNN*i):(G2obsNN*i+(G2obsNN-1))][]);
			G2obsYensQ=G2obsYensQ+G2obsNN*(((G2obsbigO-G2obsbigE).^2)./((G2obsbigE).*(1-G2obsbigE)));
		}

			// Chisquare
		decl G2obsP1out=TwoPL(repaG2,repbG2,grid);
		decl G2obsP0out=1-G2obsP1out;
		decl G2obsPout=G2obsP0out~G2obsP1out;
		decl G2obsPP0out=grp2N*sumc(G2obsP0out.*(pgrid*ones(1,grp2J)))/10;
		decl G2obsPP1out=grp2N*sumc(G2obsP1out.*(pgrid*ones(1,grp2J)))/10;
		decl G2obsExp0=G2obsPP0out';
		decl G2obsExp1=G2obsPP1out';
		decl G2obsObs0=sumc((1-grp2x))';
		decl G2obsObs1=sumc(grp2x)';
		decl G2obsChisq0=((G2obsObs0-G2obsExp0).^2)./G2obsExp0;
		decl G2obsChisq1=((G2obsObs1-G2obsExp1).^2)./G2obsExp1;
		decl G2obsChisq=G2obsChisq0+G2obsChisq1;
		
		// Replicated Discrepency Statistics

			// Observed Score Distribution
		decl G2repNC=sumc(repXG2);
		decl G2repENC=sumc(repPG2);
		decl G2repOSD=((G2repNC-G2repENC).^2)./G2repENC;
		for(decl k=0;k<grp2J;k++){
			OSDsampleG2[k][counter-1]=G2repOSD[k];
		}
			// Outfit
		decl G2repoutfitE = repPG2;
		decl G2repoutfitp0 = 1-G2repoutfitE;
		decl G2repoutfitp1 = G2repoutfitE;
		decl G2repoutfitVar = (((0 - G2repoutfitE).^2).*G2repoutfitp0) + (((1 - G2repoutfitE).^2).*G2repoutfitp1);
		decl G2repoutfitZ = (repXG2 - G2repoutfitE)./sqrt(G2repoutfitVar);
		decl G2repoutfitsq = G2repoutfitZ.^2;
		decl G2repoutfit = sumc(G2repoutfitsq)/grp2N;
		for(decl k=0;k<grp2J;k++){
			OUTsampleG2[k][counter-1]=G2repoutfit[k];
		}
			// Infit
		decl G2repinfit = sumc(G2repoutfitVar.*G2repoutfitsq)./sumc(G2repoutfitVar);
		for(decl k=0;k<grp2J;k++){
			INsampleG2[k][counter-1]=G2repinfit[k];
		}
			// Yen's Q1
		decl G2reptempmatrix=repXG2~repthG2;
		decl G2repordermatrix=sortbyc(G2reptempmatrix,grp2J);
		decl G2reporderX=dropc(G2repordermatrix,grp2J);
		decl G2reporderth=G2repordermatrix[][grp2J];
		decl G2reporderP=TwoPL(repaG2,repbG2,G2reporderth);
		decl G2repNN=grp2N/10;
		decl G2repYensQ=zeros(1,grp2J);
		for(decl i=0;i<10;i++){
			decl G2repbigO=(1/G2repNN)*sumc(G2reporderX[(G2repNN*i):(G2repNN*i+(G2repNN-1))][]);
			decl G2repbigE=(1/G2repNN)*sumc(G2reporderP[(G2repNN*i):(G2repNN*i+(G2repNN-1))][]);
			G2repYensQ=G2repYensQ+G2repNN*(((G2repbigO-G2repbigE).^2)./((G2repbigE).*(1-G2repbigE)));
		}
		for(decl k=0;k<grp2J;k++){
			Q1sampleG2[k][counter-1]=G2repYensQ[k];
		}

			// Chisquare
		decl G2repP1out=TwoPL(repaG2,repbG2,grid);
		decl G2repP0out=1-G2repP1out;
		decl G2repPout=G2repP0out~G2repP1out;
		decl G2repPP0out=grp2N*sumc(G2repP0out.*(pgrid*ones(1,grp2J)))/10;
		decl G2repPP1out=grp2N*sumc(G2repP1out.*(pgrid*ones(1,grp2J)))/10;
		decl G2repExp0=G2repPP0out';
		decl G2repExp1=G2repPP1out';
		decl G2repObs0=sumc((1-repXG2))';
		decl G2repObs1=sumc(repXG2)';
		decl G2repChisq0=((G2repObs0-G2repExp0).^2)./G2repExp0;
		decl G2repChisq1=((G2repObs1-G2repExp1).^2)./G2repExp1;
		decl G2repChisq=G2repChisq0+G2repChisq1;
		for(decl k=0;k<grp2J;k++){
			ChisqsampleG2[k][counter-1]=G2repChisq[k];
		}


		// Compute the predictive posterior p-value (ppp)
		decl G2OSDind=vecindex(G2repOSD.>G2obsOSD);
		decl G2outfitind=vecindex(G2repoutfit.>G2obsoutfit);
		decl G2infitind=vecindex(G2repinfit.>G2obsinfit);
		decl G2YensQind=vecindex(G2repYensQ.>G2obsYensQ);
		decl G2Chisqind=vecindex(G2repChisq.>G2obsChisq);
		OSDpvalG2[G2OSDind]=OSDpvalG2[G2OSDind]+1;
		OUTpvalG2[G2outfitind]=OUTpvalG2[G2outfitind]+1;
		INpvalG2[G2infitind]=INpvalG2[G2infitind]+1;
		Q1pvalG2[G2YensQind]=Q1pvalG2[G2YensQind]+1;
		ChisqpvalG2[G2Chisqind]=ChisqpvalG2[G2Chisqind]+1;

		counter=counter+1;
	}
}

print("\n","------- Posterior Predictive P-value (Sinharay & Johnson, 2003) ------");
print("\n","-------  .05 < Goodfit < .95  ------","\n");
print("\n","-------   By Item (Group 1)   ------");
print( "%c", {"Item","Infit PPP","Outfit PPP","OSD PPP","Q1 PPP","Chisq PPP"}, "%cf",{"%10.0f","%15.4f","%15.4f","%15.4f","%15.4f","%15.4f"}, range(1, grp1J, 1)'~(INpvalG1~OUTpvalG1~OSDpvalG1~Q1pvalG1~ChisqpvalG1)/chains/(it-burn));

print("\n","-------   By Item (Group 2)   ------");
print( "%c", {"Item","Infit PPP","Outfit PPP","OSD PPP","Q1 PPP","Chisq PPP"}, "%cf",{"%10.0f","%15.4f","%15.4f","%15.4f","%15.4f","%15.4f"}, range(1, grp2J, 1)'~(INpvalG2~OUTpvalG2~OSDpvalG2~Q1pvalG2~ChisqpvalG2)/chains/(it-burn));



//--- Output results to file, if option selected;
if(outputtheta==1){
	format(252);
	decl outres=sprint(resultfile);
	decl fileres=fopen(outres,"w");
	fprintln(fileres,"Parameter Means:");
	for(decl i=0;i<totalJ;i++){
		for(decl p=0;p<2;p++){
			fprint(fileres,means[i][p],"\t");
		}
		fprintln(fileres," ");
	}
	fprintln(fileres,"Parameter SDs:");
	for(decl i=0;i<totalJ;i++){
		for(decl p=0;p<2;p++){
			fprint(fileres,SDs[i][p],"\t");
		}
		fprintln(fileres," ");
	}
	fclose(fileres);
}


//--- Output thetas to file, if option selected;
if(outputtheta==1){
	format(252);
	decl outth=sprint(thetafile);
	decl fileth=fopen(outth,"w");
	fprint(fileth,"%#M",sthsum/((i-burn)*chains));
	fclose(fileth);
}

//
//--- Output thetas to file, if option selected;
decl outmc, filemc;
if(outputsamples==1){
	format(252);
	outmc=sprint(outputINsamplesG1);
	filemc=fopen(outmc,"w");
	fprint(filemc,"%#M",INsampleG1');
	fclose(filemc);

	format(252);
	outmc=sprint(outputINsamplesG2);
	filemc=fopen(outmc,"w");
	fprint(filemc,"%#M",INsampleG2');
	fclose(filemc);

	format(252);
	outmc=sprint(outputOUTsamplesG1);
	filemc=fopen(outmc,"w");
	fprint(filemc,"%#M",OUTsampleG1');
	fclose(filemc);

	format(252);
	outmc=sprint(outputOUTsamplesG2);
	filemc=fopen(outmc,"w");
	fprint(filemc,"%#M",OUTsampleG2');
	fclose(filemc);

	format(252);
	outmc=sprint(outputOSDsamplesG1);
	filemc=fopen(outmc,"w");
	fprint(filemc,"%#M",OSDsampleG1');
	fclose(filemc);

	format(252);
	outmc=sprint(outputOSDsamplesG2);
	filemc=fopen(outmc,"w");
	fprint(filemc,"%#M",OSDsampleG2');
	fclose(filemc);

	format(252);
	outmc=sprint(outputQ1samplesG1);
	filemc=fopen(outmc,"w");
	fprint(filemc,"%#M",Q1sampleG1');
	fclose(filemc);

	format(252);
	outmc=sprint(outputQ1samplesG2);
	filemc=fopen(outmc,"w");
	fprint(filemc,"%#M",Q1sampleG2');
	fclose(filemc);
}


}



