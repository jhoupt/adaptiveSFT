data { 
   int<lower=1> N;
   real intensity[N];
   int<lower=0,upper=1> correct[N];
   real<lower=0> minRT;
   real<lower=0> rt[N];
}
transformed data {
   real square_intensity[N];
   square_intensity = square(intensity);
}
parameters{
   real alpha;
   real alpha2;
   real mu;
   real<lower=0> varZ;
   real<lower=0,upper=minRT> psi;
}
transformed parameters {
   real z[2,N];

   for (tr in 1:N) { 
      z[1,tr] = mu - alpha * intensity[tr] - alpha2 * square_intensity[tr];
      z[2,tr] = mu + alpha * intensity[tr] + alpha2 * square_intensity[tr];
   }


}
model {
   varZ ~ inv_gamma(1,.1);
   mu ~ normal(0,1);
   alpha ~ normal(0,2);
   alpha2 ~ normal(0,1);

   // psi has improper flat prior on positive reals
   for ( tr in 1:N) { 
      if ( correct[tr] ) {
         target += lognormal_lpdf(rt[tr] - psi | z[1,tr], varZ);
         target += lognormal_lccdf(rt[tr] - psi | z[2,tr], varZ);
      } 
      else { 
         target += lognormal_lpdf(rt[tr] - psi | z[2,tr], varZ);
         target += lognormal_lccdf(rt[tr] - psi | z[1,tr], varZ);
      }
   }
   
   
}
