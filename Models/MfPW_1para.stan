functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
    
    real x00 = y[1];
    real x10 = y[2];
    real x11 = y[3];
    real x01 = y[4];
    real z00 = y[5];
    
    real y1[5];
    real dx00_dt = y1[1];
    real dx10_dt = y1[2];
    real dx11_dt = y1[3];
    real dx01_dt = y1[4];
    real dz00_dt = y1[5];
    
    real lamda = theta[1];
    real tau = theta[2];
    real sigma = theta[3];
    real beta = theta[4];
    real k1 = theta[5];

    if((0<t)|| (0==t) && (t<1)){
      lamda=k1;}
    
    if((1<t) || (1==t) && (t<2)){
      lamda=k1;}
    
    if((2<t) || (2==t) && (t<3)){
      lamda=k1;}
    
    if((3<t) || (3==t) && (t<4)){
      lamda=k1;}
    
    if((4<t) || (4==t) && (t<5)){
      lamda=k1;}
    
    if((5<t) || (5==t) && (t<6)){
      lamda=k1;}
    
    if((6<t) || (6==t) && (t<7)){
      lamda=k1;}
    
    if((7<t) || (7==t) && (t<8)){
      lamda=k1;}
    
    if((8<t) || (8==t) && (t<9)){
      lamda=k1;}
    
    if((9<t) || (9==t) && (t<10)){
      lamda=k1;}
    
    if((10<t) || (10==t) && (t<11)){
      lamda=k1;}
    
    if((11<t) || (11==t) && (t<12)){
      lamda=k1;} 
    
    if((12<t) || (12==t) && (t<13)){
      lamda=k1;}
    
    if((13<t)|| (13==t) && (t<14)){
      lamda=k1;}
    
    if((14<t) || (14==t) && (t<15)){
      lamda=k1;}
    
    if((15<t) || (15==t) && (t<16)){
      lamda=k1;}
    
    if((16<t) || (16==t) && (t<17)){
      lamda=k1;}
    
    if((17<t) || (17==t) && (t<18)){
      lamda=k1;}
    
    if((18<t) || (18==t) && (t<19)){
      lamda=k1;}
    
    if((19<t) || (19==t) && (t<20)){
      lamda=k1;}
    
    if((20<t) || (20==t) && (t<21)){
      lamda=k1;}
    
    if((21<t) || (21==t) && (t<22)){
      lamda=k1;}
    
    if((22<t) || (22==t) && (t<23)){
      lamda=k1;}
    
    if((23<t) || (23==t) && (t<24)){
      lamda=k1;}
    
    if((24<t) || (24==t) && (t<25)){
      lamda=k1;} 
    
    if((25<t) || (25==t) && (t<26)){
      lamda=k1;}  
    
    if((26<t) || (26==t) && (t<27)){
      lamda=k1;}
    
    if((27<t) || (27==t) && (t<28)){
      lamda=k1;}
    
    if((28<t) || (28==t) && (t<29)){
      lamda=k1;}
    
    if((29<t) || (29==t) && (t<30)){
      lamda=k1;}
    
    if((30<t) || (30==t) && (t<31)){
      lamda=k1;}
    
    if((31<t) || (31==t) && (t<32)){
      lamda=k1;}
    
    if((32<t) || (32==t) && (t<33)){
      lamda=k1;}
    
    if((33<t) || (33==t) && (t<34)){
      lamda=k1;}
    
    if((34<t) || (34==t) && (t<35)){
      lamda=k1;}
    
    if((35<t) || (35==t) && (t<=36)){
      lamda=k1;}
    
    dx00_dt = -lamda*x00;
    dx10_dt =  lamda*x00-tau*x10;
    dx11_dt =  tau*x10 -sigma*x11;
    dx01_dt =  sigma*x11-beta*x01;
    dz00_dt =  beta * x01;
    
    return {dx00_dt, dx10_dt, dx11_dt,dx01_dt,dz00_dt};
  }
}


data {
  int<lower=1> n_week;
  int<lower=1> n_week2;
  real t0;
  real t[n_week];
  int t2[n_week2];
  int y[n_week2,4];
}
transformed data {
  int x_i[0];
  real x_r[0];
}
parameters {
  
  real <lower= 0,upper=1>  lamda;
  real <lower= 0,upper=5>  tau;
  real <lower= 0, upper=1> sigma;
  real <lower= 0, upper=1> beta;
  real <lower= 0, upper=1> y00;
  
  real <lower= 0, upper=1> k1;

  real <lower= 0, upper=1> k01;
  real <lower= 0, upper=1> k11;
  real <lower= 0, upper=1> k10;
  
  
}

transformed parameters{
  
  simplex[4] out[n_week];
  
  real<lower= 0, upper=1> y_init[5];
  real<lower= 0, upper=1> y_hat[n_week,5];
  
  {
    real theta[5];
    theta[1] = lamda;
    theta[2] = tau;
    theta[3] = sigma; 
    theta[4] = beta;
    theta[5]= k1;

    y_init[1] = y00;
    y_init[2] = k10*(1-y00);
    y_init[3] = k11*(1-y00-(k10*(1-y00)));
    y_init[4] = k01*(1-y00-(k10*(1-y00))-(k11*(1-y00-(k10*(1-y00)))));
    y_init[5] = 1-y00-k10*(1-y00)-k11*(1-y00-(k10*(1-y00)))-k01*(1-y00-(k10*(1-y00))-(k11*(1-y00-(k10*(1-y00)))));
    
    
    y_hat = integrate_ode_rk45(sir, y_init, t0, t, theta, x_r, x_i);
    
  }
  
  for(i in 1:n_week){
    out[i] = to_vector([(y_hat[i,1]+y_hat[i,5]),y_hat[i,2],y_hat[i,3],y_hat[i,4]]);
  } 
}

model {
  y00~beta(8,2); 
  
  k01~uniform(0,1);
  k11~uniform(0,1);
  k10~uniform(0,1);
  
  k1~uniform(0,1);

  tau ~ gamma(4,3);
  sigma~ uniform(0,1);
  beta ~uniform(0,1);
  
  for(i in 1:n_week2){
    target +=multinomial_lpmf(y[i,] | out[t2[i]]);
  }
  
}







