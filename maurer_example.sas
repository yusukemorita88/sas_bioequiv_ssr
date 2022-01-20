/*###########################################################################
This is the example SAS code to calculate values in Case Study (Table 6)
of the Research Article :

  Maurer, W, Jones, B, Chen, Y. 
  Controlling the type I error rate in two-stage sequential adaptive designs 
  when testing for average bioequivalence. 
  Statistics in Medicine. 2018; 37: 1587– 1607. 
  https://doi.org/10.1002/sim.7614

2021.09.20
  The R code in supplementary information (sim7614-sup-0001-supplementary.pdf)
      was ported into SAS code by Yusuke Morita;

############################################################################*/

/*#################################################################
## case study
## preliminaries
## (setting initial parameter values and calculating
## significance levels for standard and maximum combination tests)
##################################################################*/

* store required modules in IML;
%include ".\probmvt.sas"; *specify the full path for probmvt.sas;

* ABE limits (on ratio scale);
%let theta1 = 0.8;
%let theta2 = 1.25;

* assumed true ratio of Test to Reference when planning trial;
%let planned_BE_ratio = 0.95;

* significance level;
%let alpha = 0.05;

* planned power;
%let power_min = 0.8;

/*####################################################################
## significance levels for standard combination test for weight=0.25
####################################################################*/
proc iml;
    * function to solve for critical value (equal alphas in stages 1 and 2);
    start BVNCDF(x);
        w1_star = 0.25;
        f = PROBBNRM(x, x, sqrt(w1_star)) - (1 - &alpha.);
        return (f);
    finish BVNCDF;
    
    * get critical value of standard combination test;
    z_crit_stage1_star = froot("BVNCDF", {0, 5});
    
    * convert this to a significance level;
    alpha_stage1_star = 1 - probnorm(z_crit_stage1_star);
    
    print z_crit_stage1_star alpha_stage1_star;
quit;

/*###############################################################
## Calculation of critical values for maximum combination test
###############################################################*/
proc iml;
    call randseed(12345678);
    randN = 100000000/2; /* number of random variables, recommend >= 100,000,000 */
    
    * choice of weights;
    w1 = 0.5; 
    w2 = 1 - w1;
    w1_star = 0.25;
    w2_star = 1 - w1_star;
    root_w1 = sqrt(w1);
    root_w1_star = sqrt(w1_star);
    rho=sqrt(w1 * w1_star) + sqrt(w2 * w2_star);
    
    /* correlation matrix */
    Sigma = (1              || root_w1  || root_w1_star )//
            (root_w1        || 1        || rho          )//
            (root_w1_star   || rho      || 1            );
    mean = {0 0 0};
    Z = randnormal(randN, mean, Sigma);    /* sample from MVN(0, Sigma) */
    create MVNorm3 from Z[c={Z1 Z0a Z0b}];  append from Z;  close;

    /*#########################################################
    ## function used in froot to solve for critical value
    ## of maximum combination test
    #########################################################*/
    start TVNCDF(x);
        use MVNorm3;
        read all into TRIVN;
        ind = (TRIVN[,1] < x) & (TRIVN[,2] < x) & (TRIVN[,3] < x);
        prob = ind[:,];
        f = (prob) - (1 - &alpha.);
        return (f);
    finish TVNCDF;
    
    * get critical value of robust combination test for z_max1 and z_max2;
    z_crit_max = froot("TVNCDF", {0, 5});
    
    * convert to nominal significance level;
    alpha_max = 1 - probnorm(z_crit_max);
    
    * store to macro variable;
    call symput("z_crit_max", char(z_crit_max));
    call symput("alpha_max", char(alpha_max));
    print z_crit_max alpha_max;
quit;

/*######################
## single stage design
######################*/
* sample size calculation for CV=0.3, alpha=0.05 and target power=0.8;
%let intrasubjectCV = 0.3;*macro vaiable for the intra-subject CV;
%let sigma2 = %sysfunc(log(1+&intrasubjectCV.**2));* macro variable for the within-subject variance;
%let std_derived = %sysevalf((&sigma2./2)**0.5); * macro variable for the derived common standard deviation;
%let log_pt_8 = %sysfunc(log(&theta1.)); * macro variable for log(0.8);
%let log_1_pt_25 = %sysfunc(log(&theta2.)); * macro variable for log(1.25);
%let log_true_gmr = %sysfunc(log(&planned_BE_ratio.)); * macro variable for log(planned_BE_ratio);

proc power;
     twosamplemeans test=equiv_diff alpha= &alpha.
     lower=&log_pt_8. upper=&log_1_pt_25. std=&std_derived.
     meandiff=&log_true_gmr.
     ntotal=.
     power = &power_min.;
run; 

* sample size calculation for CV=0.4, alpha=0.05 and target power=0.8;
%let intrasubjectCV = 0.4;*macro vaiable for the intra-subject CV;
%let sigma2 = %sysfunc(log(1+&intrasubjectCV.**2));* macro variable for the within-subject variance;
%let std_derived = %sysevalf((&sigma2./2)**0.5); * macro variable for the derived common standard deviation;

proc power;
     twosamplemeans test=equiv_diff alpha= &alpha.
     lower=&log_pt_8. upper=&log_1_pt_25. std=&std_derived.
     meandiff=&log_true_gmr.
     ntotal=.
     power = &power_min.;
run; 

/*#################################################################
## Simulate data (within-subject differences) for first stage,
## estimate model parameters,
## apply TOST procedure and calculate achieved power of first stage
##################################################################*/
*from Table 3;
data stage1;
    input GRP $ DIFF;
cards;
RT -0.1105 
RT 0.6398 
RT 0.4655 
RT -1.2661 
RT -0.5989 
RT -0.4064 
RT -0.4232 
RT 0.4906 
RT 0.1440 
RT -0.2435 
TR -0.5600
TR -0.2480
TR -0.0586
TR -0.1810
TR 0.9407
TR 0.2048 
TR 0.0457 
TR -0.1744 
TR -0.3595
TR -0.0694
;
run;

data stage1a;
    set stage1;
    USUBJID = _N_;
run;

proc univariate data = stage1a;
    class GRP;
    var DIFF;
    output out = res n = N mean = MEAN var = VAR;
run;

data res1;
    merge
        res(where=(GRP="TR") rename=(N=N_TR MEAN=MEAN_TR VAR=VAR_TR))
        res(where=(GRP="RT") rename=(N=N_RT MEAN=MEAN_RT VAR=VAR_RT))
    ;
run;

data interim;
    set res1;
    
    *size of first stage;
    n_stage1 = N_TR + N_RT;
    
    * estimate of mean difference;
    delta_hat_stage1 = ( MEAN_TR - MEAN_RT ) / 2;
    
    * error degrees of freedom;
    errdf_stage1 = n_stage1 - 2;
    
    * unblinded estimate of variance;
    var_hat_stage1 = 0.5 * (VAR_TR*(N_TR-1) + VAR_RT*(N_RT-1)) /errdf_stage1;
    
    * estimated standard deviation;
    sigma_hat_stage1 = sqrt(var_hat_stage1);
    
    * estimated CV;
    cv_hat_stage1 = sqrt(exp(var_hat_stage1)-1);
    
    * standard error of treatment difference on log scale;
    stderr_stage1 = sqrt(2*var_hat_stage1/n_stage1);
    
    * 90percent confidence interval;
    t_crit = tinv(0.95, errdf_stage1);*t-critical value;
    ci_low = delta_hat_stage1 - t_crit * stderr_stage1;
    ci_upp = delta_hat_stage1 + t_crit * stderr_stage1;
    
    ci_low_exp = exp(ci_low);
    ci_upp_exp = exp(ci_upp);
    
    * test statistics for stage 1;
    t1_stage1 = (delta_hat_stage1 - log(0.8)) / stderr_stage1;
    t2_stage1 = (log(1.25) - delta_hat_stage1) / stderr_stage1;
    
    * p-values for stage 1 (based on t-tests);
    p_val_t1_stage1 = 1 - probt(t1_stage1, errdf_stage1);
    p_val_t2_stage1 = 1 - probt(t2_stage1, errdf_stage1);
    
    * convert p-values into z-statistics;
    z1_stage1 = probit(1 - p_val_t1_stage1);
    z2_stage1 = probit(1 - p_val_t2_stage1);
    
    * apply futility rule;
    futile = 0;
    if(exp(ci_low) >= 0.95 and exp(ci_low) <= 1.05) then futile = 1;
    if(exp(ci_upp) >= 0.95 and exp(ci_upp) <= 1.05) then futile = 1;
run;

/*#############################################
## TOST for first stage using z-statistics and
## z critical value based on alpha_stage1
#############################################*/
data tost1st;
    set interim;
    
    z_crit_stage1 = 1.9374;%*&z_crit_max.;
    * can null hypothesis be rejected at left side?;
    z_nrejlow_stage1 = (z1_stage1 >= z_crit_stage1); *1:rejected, 0:accepted;
    * can null hypothesis be rejected at right side;
    z_nrejupp_stage1 = (z2_stage1 >= z_crit_stage1); *1:rejected, 0:accepted;

    * can reject null hypothesis at both sides?;
    * i.e., ABE accepted;
    z_nABE_passed_stage1 = (z_nrejlow_stage1 & z_nrejupp_stage1 & futile = 0); *1:ABE accepted, 0:ABE rejected;
    
    keep z: futile;
run;

/*#######################################################
## sample size calculation for Method B using CV=0.3682,
## alpha=0.0294 and delta=log(0.95)
########################################################*/
%let intrasubjectCV = 0.36821449;*macro vaiable for the intra-subject CV (cv_hat_stage1);
%let sigma2 = %sysfunc(log(1+&intrasubjectCV.**2));* macro variable for the within-subject variance;
%let std_derived = %sysevalf((&sigma2./2)**0.5); * macro variable for the derived common standard deviation;
%let log_pt_8 = %sysfunc(log(&theta1.)); * macro variable for log(0.8);
%let log_1_pt_25 = %sysfunc(log(&theta2.)); * macro variable for log(1.25);
%let log_true_gmr = %sysfunc(log(&planned_BE_ratio.)); * macro variable for log(planned_BE_ratio);

proc power;
     twosamplemeans test=equiv_diff alpha= 0.0294
     lower=&log_pt_8. upper=&log_1_pt_25. std=&std_derived.
     meandiff=&log_true_gmr.
     ntotal=.
     power = &power_min.;
run; 

/*####################################################
## achieved power of first stage for
## maximum combination test with weights = (0.5,0.25)
####################################################*/
%let n_stage1=20;

ods output output=achieved_power_stage1;
proc power;
     twosamplemeans test=equiv_diff alpha= 0.0264 %*&alpha_max.;
     lower=&log_pt_8. upper=&log_1_pt_25. std=&std_derived.
     meandiff=&log_true_gmr.
     ntotal= &n_stage1.
     power = .;
run;

/*####################################################
## conditional errors for the maximum combination test
####################################################*/
data conderr;
    set interim;
    * nominal significance levels critical values for maximum combination test;
    z_crit_max = 1.9374;%*&z_crit_max.;
    
    w1 = 0.5;
    w2 = 1 - w1;
    w1_star = 0.25;
    w2_star = 1 - w1_star;

    
    * calculate value within Phi function for H{01};
    val_11 = (z_crit_max - sqrt(w1) * z1_stage1) / sqrt(w2);
    val_12 = (z_crit_max - sqrt(w1_star) * z1_stage1) / sqrt(w2_star);
    val_min1 = min(val_11, val_12);
    alpha_c1_max = 1 - probnorm(val_min1);
    
    * calculate values within Phi function for H{02};
    val_21 = (z_crit_max - sqrt(w1) * z2_stage1) / sqrt(w2);
    val_22 = (z_crit_max - sqrt(w1_star) * z2_stage1) / sqrt(w2_star);
    val_min2 = min(val_21, val_22);
    alpha_c2_max = 1 - probnorm(val_min2);
    
    t_crit1 = tinv((1 - alpha_c1_max), errdf_stage1);
    t_crit2 = tinv((1 - alpha_c2_max), errdf_stage1);
    
    keep val_: alpha_: t_:;
run;

/*####################################################
## planned Type II error and target conditional power
####################################################*/
data condpower;
    set achieved_power_stage1;
    planned_beta = 1 - &power_min.;
    beta1_hat = 1 - power; * power = power_stage1_alpha_max;
    cond_power = (beta1_hat - planned_beta) / beta1_hat;
    keep planned_beta beta1_hat power cond_power;
    format power best12.;
run;

/*###################
## adaptive planning
###################*/
data adaptplan;
    set interim;
    adapt_planned_BE_ratio = &planned_BE_ratio.;
    adapt_planned_delta = log(adapt_planned_BE_ratio);
    delta_planned_ratio = log(&planned_BE_ratio.);
    if(sign(adapt_planned_delta) ^= sign(delta_hat_stage1)) then do;
        adapt_planned_BE_ratio = 1 / &planned_BE_ratio.;
    end;
    keep adapt_planned: delta_:;
run;


/*####################################################################
## calculate sample size of second stage to achieve conditional power
####################################################################*/
*copy variable to macro variable;
proc sql noprint;
    select adapt_planned_BE_ratio
    into  :adapt_planned_BE_ratio
    from adaptplan;
    
    select alpha_c1_max, alpha_c2_max 
    into  :alpha_c1_max,:alpha_c2_max
    from conderr;
    
    select n_stage1, var_hat_stage1
    into  :n_stage1,:var_hat_stage1
    from interim;

    select cond_power
    into  :cond_power
    from condpower;
quit;

/*################################################
## power function for unequal significance levels
################################################*/
proc iml;
    load ; *load modules in probmvt.sas;

    lnGmrPlanned = log(&adapt_planned_BE_ratio.);
    ln1p25 = log(&theta2.);
    variance1 = &var_hat_stage1.;
    
    COVAR = {1 -1, 
            -1  1};
    acheived = 0;
    
    do N = %eval(&n_stage1.+2) to 200 by 2 until(acheived);

        DIM = 2;
        DF = N - 2;
        
        nc1 = (ln1p25 + lnGmrPlanned) /sqrt(variance1 * 2 / N);
        nc2 = (ln1p25 - lnGmrPlanned) /sqrt(variance1 * 2 / N); 
        DELTA = j(1, 2);
        DELTA[1] = nc1;
        DELTA[2] = nc2;
    
        INFIN = J(1, DIM, 1);
        tc1 = tinv(1 - &alpha_c1_max., DF);
        tc2 = tinv(1 - &alpha_c2_max., DF);
        LOWER = j(1, 2);
        LOWER[1] = tc1;
        LOWER[2] = tc2;     
        UPPER = j(1, 2, 999);*infinity:999 is ignored;
        
        MAXPTS = 2000*DIM*DIM*DIM;
        ABSEPS = .0000001;
        RELEPS = 0;
  
        RUN MVN_DIST( DIM, DF, DELTA, LOWER, UPPER, INFIN, COVAR, MAXPTS, ABSEPS, RELEPS,  ERROR, condpower, NEVALS, INFORM );

        if (condpower >= &cond_power.) then acheived = 1;
        print N condpower;
        if acheived = 1 then do;
            create SSR from N[c={N}];  
            append from N;  
            close;
        end;
    end;
    print N;
quit;

/*####################################################
## simulated data (as in paper) for second stage
####################################################*/
*from Table 5;
data stage2;
    input GRP $ DIFF;
cards;
RT 0.0322 
RT -0.0503 
RT 1.0717 
RT 0.5169 
RT -0.0844 
RT 0.2189 
RT -0.7633 
RT 0.2539 
RT 0.0063 
RT -0.1603 
RT -0.4497 
RT -0.1348 
RT 0.5409 
RT -0.4727 
RT 0.0973 
RT 0.0401 
RT 0.0482
RT -0.4781
TR -0.7412
TR 1.0219
TR 0.0050
TR 0.2836
TR -0.8620
TR -0.2195
TR -0.6296
TR -0.4350
TR -0.0490
TR 0.4293
TR -0.1376
TR 0.1901
TR -0.5432
TR 1.0669
TR 0.6068
TR -0.4433
TR -0.0512
TR 0.2593
;
run;


data stage2a;
    set stage2;
    USUBJID = _N_;
run;

proc univariate data = stage2a;
    class GRP;
    var DIFF;
    output out = res2 n = N mean = MEAN var = VAR;
run;

data res2a;
    merge
        res2(where=(GRP="TR") rename=(N=N_TR MEAN=MEAN_TR VAR=VAR_TR))
        res2(where=(GRP="RT") rename=(N=N_RT MEAN=MEAN_RT VAR=VAR_RT))
    ;
run;

/*#######################################
## estimate parameters for second stage
########################################*/
data res2nd;
    set res2a;
    * number of data values for second stage;
    n_stage2 = N_TR + N_RT;
    
    * estimate of mean difference for second stage;
    delta_hat_stage2 = ( MEAN_TR - MEAN_RT ) / 2;
    
    * df for second stage;
    errdf_stage2 = n_stage2 - 2;
    
    * unblinded estimate of variance;
    var_hat_stage2 = 0.5 * (VAR_TR*(N_TR-1) + VAR_RT*(N_RT-1)) /errdf_stage2;

    * estimated CV for second stage;
    cv_hat_stage2 = sqrt(exp(var_hat_stage2)-1);
    
    * estimated standard deviation;
    sigma_hat_stage2 = sqrt(var_hat_stage2);
    
    * standard error of treatment difference for stage 2 only;
    stderr_stage2 = sqrt(2*var_hat_stage2/n_stage2);
    
    
    /*#######################
    ## TOST for second stage
    ########################*/
    * one-sided t-test values for stage 2;
    t1_stage2 = (delta_hat_stage2 - log(&theta1.)) / stderr_stage2;
    t2_stage2 = (log(&theta2.) - delta_hat_stage2) / stderr_stage2;
    
    * p-values for stage 2 (based on t-tests);
    p_val_t1_stage2 = 1 - probt(t1_stage2, errdf_stage2);
    p_val_t2_stage2 = 1 - probt(t2_stage2, errdf_stage2);
    
    * convert these into z-statistics for use in the maximum combination test;
    z1_stage2 = probit(1 - p_val_t1_stage2);
    z2_stage2 = probit(1 - p_val_t2_stage2);
    
run;

/*#########################################################
## TOST at end of trial using the maximum combination test
## using weights (0.5,0.5) and (0.25,0.75)
#########################################################*/
data final;
    merge interim res2nd;
    w = 0.5;
    w_star = 0.25;
    z_crit_max = 1.9374;%*&z_crit_max.;
    
    * z-statistics;
    z1_w = sqrt(w) * z1_stage1 + sqrt(1 - w) * z1_stage2;
    z2_w = sqrt(w) * z2_stage1 + sqrt(1 - w) * z2_stage2;
    
    z1_w_star = sqrt(w_star) * z1_stage1 + sqrt(1 - w_star) * z1_stage2;
    z2_w_star = sqrt(w_star) * z2_stage1 + sqrt(1 - w_star) * z2_stage2;
    
    z1_max = max(z1_w, z1_w_star);
    z2_max = max(z2_w, z2_w_star);
    
    * reject null hypothesis at left side?;
    nrejlow = (z1_max >= z_crit_max);
    
    * reject null hypothesis at right side?;
    nrejupp = (z2_max >= z_crit_max);
    
    * reject null hypothesis at both sides?,;
    * i.e., is ABE accepted?;
    z_passed = (nrejlow & nrejupp); *1:ABE accepted, 0:ABE rejected;

    keep w w_star z_crit_max z1_w z1_w_star z1_max z2_w z2_w_star z2_max nrejlow nrejupp z_passed;
run;

*EOF;