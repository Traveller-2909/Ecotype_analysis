//Lifelong growth dynamics of pike ecotypes
//Hierarchical design with three levels:
//   (1) multiple data (otolith annual growth marks) per fish 
//   (2) several fish per ecotype
//   (3) several ecotypes
// Main objective: Looking for differences in Linf and k between the ecotypes
functions {
	real gamm_a(real phi_inv) {
		return phi_inv;
	}
	real gamm_b(real mu, real phi_inv) {
		return phi_inv / mu;
	}
	vector gamm_a_v(real phi_inv, int n) {
		return rep_vector(phi_inv, n);
	}
	vector gamm_b_v(vector mu, real phi_inv) {
		return phi_inv ./ mu;
	}
}
data {
	int <lower=0> N; //Number of observations
	
	//data vectors
	vector <lower=0, upper=13> [N] Agei; //Age in years
	vector <lower=0> [N] Radmm; //otolith growth increments (in mm)
	
	//groups
	int<lower=1> n_fish; //number of individual fish
	int<lower=1> n_eco; //number of ecotypes
	int<lower=1, upper=n_fish> fish_id [N];
	int <lower = 1, upper = n_eco> eco_id_fish [n_fish];
}
parameters {

	// fish-level
	vector <lower=0> [n_fish] LInf_i;  //Mean maximum adult size at fish level
	vector <lower=0> [n_fish] tZero_i; //Age at wich size = 0 at fish level (size can´t be zero, so has to be negative)
	vector <lower=0> [n_fish] k_i;     //Growth completion parameter at fish level

	// ecotype level
	vector <lower=0> [n_eco] LInf_eco; //Mean maximum adult size at ecotype level
	vector <lower=0> [n_eco] tZero_eco;//Age at wich size = 0 at ecotype level
	vector <lower=0> [n_eco] k_eco;    //Growth completion parameter at ecotype level

	//hyperprameters
	real <lower=0> LInf;     //Baseline expectation for mean maximum adult size
	real <lower=0> tZero;    //Baseline expectation for age at wich size = 0
	real <lower=0> k;        //Baseline expectation for growth completion parameter

	// dispersion parameters
	real <lower=0> phi;
	real <lower=0> LInf_phi; 
	real <lower=0> k_phi;    
	real <lower=0> tZero_phi;
	vector <lower=0> [n_eco] LInf_eco_phi; 
	vector <lower=0> [n_eco] k_eco_phi;    
	vector <lower=0> [n_eco] tZero_eco_phi;
}
transformed parameters {
	vector[N] mu;

   // VBGM Expectation for Length at age (fish level)
	for(i in 1:N) {
		int fid = fish_id[i];
		mu[i] = LInf_i[fid] * (1 - exp(-k_i[fid] * (Agei[i] + tZero_i[fid])));
	}
}
model {
	// observation level
	Radmm ~ gamma(gamm_a_v(phi, N), gamm_b_v(mu, phi));

	// fish-level
	for(i in 1:n_fish) {
		int eid = eco_id_fish[i];
		LInf_i[i] ~ gamma(gamm_a(LInf_eco_phi[eid]), gamm_b(LInf_eco[eid], LInf_eco_phi[eid]));
		tZero_i[i] ~ gamma(gamm_a(tZero_eco_phi[eid]), gamm_b(tZero_eco[eid], tZero_eco_phi[eid]));
		k_i[i] ~ gamma(gamm_a(k_eco_phi[eid]), gamm_b(k_eco[eid], k_eco_phi[eid]));
	}

	// ecotype-level
	LInf_eco ~ gamma(gamm_a(LInf_phi), gamm_b(LInf, LInf_phi));
	k_eco ~ gamma(gamm_a(k_phi), gamm_b(k, k_phi));
	tZero_eco ~ gamma(gamm_a(tZero_phi), gamm_b(tZero, tZero_phi));
	
	//hyperpriors
	// LInf ~ normal(2.3, 0.5); //2.3 mm was the largest otolith radius in sample
	// k ~ normal(0.1, 0.05);
	// tZero ~ normal(1, 0.3);

	// dispersion parameters
	LInf_phi ~ cauchy(0, 200);
	tZero_phi ~ cauchy(0, 200);
	k_phi ~ cauchy(0, 200);
	LInf_eco_phi ~ cauchy(0, 200);
	k_eco_phi ~ cauchy(0, 200);
	tZero_eco_phi ~ cauchy(0, 200);
	phi ~ cauchy(0, 200);

	// proposed revisions
	LInf ~ normal(0, 2.5);
	k ~ normal(0, 0.5);
	tZero ~ normal(0, 1);
}
