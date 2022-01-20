/*###########################################################################
This is the example SAS code of the Maximum Combination Test 
     proposed in the Research Article :

  Maurer, W, Jones, B, Chen, Y. 
  Controlling the type I error rate in two-stage sequential adaptive designs 
  when testing for average bioequivalence. 
  Statistics in Medicine. 2018; 37: 1587– 1607. 
  https://doi.org/10.1002/sim.7614

2021.10.30: Yusuke Morita, new for SAS User Group Meeting in Japan 2021;

############################################################################*/
**********************************************************************;
* configulation;
**********************************************************************;

* ABE limits (on ratio scale);
%let theta1 = 0.8;
%let theta2 = 1.25;

* assumed true ratio of Test to Reference when planning trial;
%let planned_BE_ratio = 0.95;

* significance level(one-sided);
%let alpha = 0.05;

* planned power;
%let power_min = 0.8;

* weights for maximum combination test;
%let w = 0.5;
%let w_star = 0.25;

**********************************************************************;
* Calculation of critical values for maximum combination test;
**********************************************************************;
proc iml;
    call randseed(12345678);
    randN = 100000000/2; /* number of random variables, recommend >= 100,000,000 */
    
    * choice of weights;
    w1 = &w.; 
    w2 = 1 - w1;
    w1_star = &w_star.;
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

**********************************************************************;
* Sample Data based on Maurer's Table 3 and Table 5;
**********************************************************************;
data maurer_table3_table5;
    input STAGE TRTSEQA $ DIFF;
cards;
1 RT -0.1105 
1 RT 0.6398 
1 RT 0.4655 
1 RT -1.2661 
1 RT -0.5989 
1 RT -0.4064 
1 RT -0.4232 
1 RT 0.4906 
1 RT 0.1440 
1 RT -0.2435 
1 TR -0.5600
1 TR -0.2480
1 TR -0.0586
1 TR -0.1810
1 TR 0.9407
1 TR 0.2048 
1 TR 0.0457 
1 TR -0.1744 
1 TR -0.3595
1 TR -0.0694
2 RT 0.0322 
2 RT -0.0503 
2 RT 1.0717 
2 RT 0.5169 
2 RT -0.0844 
2 RT 0.2189 
2 RT -0.7633 
2 RT 0.2539 
2 RT 0.0063 
2 RT -0.1603 
2 RT -0.4497 
2 RT -0.1348 
2 RT 0.5409 
2 RT -0.4727 
2 RT 0.0973 
2 RT 0.0401 
2 RT 0.0482
2 RT -0.4781
2 TR -0.7412
2 TR 1.0219
2 TR 0.0050
2 TR 0.2836
2 TR -0.8620
2 TR -0.2195
2 TR -0.6296
2 TR -0.4350
2 TR -0.0490
2 TR 0.4293
2 TR -0.1376
2 TR 0.1901
2 TR -0.5432
2 TR 1.0669
2 TR 0.6068
2 TR -0.4433
2 TR -0.0512
2 TR 0.2593
;
run;

data stage1plus2;
    call streaminit(12345678);
    set maurer_table3_table5;
    USUBJID = _n_;
    _rand = rand('Uniform');
    TRTA1 = substr(TRTSEQA, 1, 1);
    TRTA2 = substr(TRTSEQA, 2, 1);
    lnCmax1 = _rand + diff; * diff = period1 - period2;
    lnCmax2 = _rand ;
    keep USUBJID STAGE TRTSEQA TRTA: lnCmax:;
run;

**********************************************************************;
* Interim Analysis;
**********************************************************************;
ods output Statistics=stat1 ConfLimits=cl1 EquivLimits=ev1 EquivTests=et1; 
proc ttest data=stage1plus2 tost(%sysfunc(log(&theta1.)), %sysfunc(log(&theta2.))) alpha=&alpha. order=internal;
    where STAGE = 1;
    var lnCmax1 lnCmax2 / crossover=(TRTA1 TRTA2) ;
run;

*---------------------------------------;
* check futility;
*---------------------------------------;
data cl90pct1;
    set ev1;
    where Treatment = "Diff (1-2)" and Method="Pooled";
    delta_hat_stage1 = exp(-mean); * minus to convert diff(R - T) to diff(T - R);
    lower90_stage1 = exp(-upperCLMean);
    upper90_stage1 = exp(-lowerCLMean);
run;

*---------------------------------------;
* apply TOST at end of Stage 1;
*---------------------------------------;
data et1a;
    set et1;
    where Treatment = "Diff (1-2)" and Method="Pooled" and Test ^= "Overall";
    if Test = "上限" then Test = "Upper";
    if Test = "下限" then Test = "Lower";
run;

proc transpose data = et1a out = et1b prefix=P1st_;
    by DF; * to keep in dataset;
    var probt;
    id Test;
run;

data tost1st;
    set et1b;
    * convert p-values into z-statistics;
    z1_stage1 = probit(1 - P1st_Lower);
    z2_stage1 = probit(1 - P1st_Upper);
    
    z_crit_stage1 = &z_crit_max.;
    
    * can null hypothesis be rejected at left side?;
    z_nrejlow_stage1 = (z1_stage1 >= z_crit_stage1); *1:rejected, 0:accepted;
    
    * can null hypothesis be rejected at right side;
    z_nrejupp_stage1 = (z2_stage1 >= z_crit_stage1); *1:rejected, 0:accepted;

    * can reject null hypothesis at both sides?;
    * i.e., ABE accepted;
    z_nABE_passed_stage1 = (z_nrejlow_stage1 & z_nrejupp_stage1); *1:ABE accepted, 0:ABE rejected;
    
    keep p1st_: z: DF;
run;

*--------------------------------------------------;
* check achieved power of first stage1 for futility;
*--------------------------------------------------;
data iscv1;
    set stat1(where=(Sequence ^= "Both" and Treatment = "Diff (1-2)")) end=last;
    retain _N _StdDev;
    
    if last then do;
        var_hat_stage1 = 0.5 * ((_N-1)*_StdDev**2 + (N-1)*StdDev**2) / (N+_N-2);
        cv_hat_stage1 = sqrt(exp(var_hat_stage1)-1);
        call symputx('intrasubjectCV', cv_hat_stage1);*macro vaiable for the intra-subject CV;
        n_stage1 = _N + N;
        call symputx('n_stage1', n_stage1);*macro vaiable for N of stage 1;
        output;
     end;
    _N = N;
    _StdDev = StdDev;
run;

%let sigma2 = %sysfunc(log(1+&intrasubjectCV.**2));* macro variable for the within-subject variance;
%let std_derived = %sysevalf((&sigma2./2)**0.5); * macro variable for the derived common standard deviation;
%let log_pt_8 = %sysfunc(log(&theta1.)); * macro variable for log(0.8);
%let log_1_pt_25 = %sysfunc(log(&theta2.)); * macro variable for log(1.25);
%let log_true_gmr = %sysfunc(log(&planned_BE_ratio.)); * macro variable for log(planned_BE_ratio);

ods output output=achieved_power_stage1;
proc power;
     twosamplemeans test=equiv_diff alpha= &alpha_max.
     lower=&log_pt_8. upper=&log_1_pt_25. std=&std_derived.
     meandiff=&log_true_gmr.
     ntotal= &n_stage1.
     power = .;
run;

*if power >= 0.80 stop (non-binding);

**********************************************************************;
* Sample Size Re-Estimation;
**********************************************************************;
*--------------------------------------------------;
* conditional error rates for the maximum combination test;
*--------------------------------------------------;
data conderr;
    set tost1st;
    * nominal significance levels critical values for maximum combination test;
    z_crit_max = &z_crit_max.;
    w1 = &w.;
    w2 = 1 - w1;
    w1_star = &w_star.;
    w2_star = 1 - w1_star;
    errdf_stage1 = DF;

    * calculate value within Phi function for H{01};
    val_11 = (z_crit_max - sqrt(w1) * z1_stage1) / sqrt(w2);
    val_12 = (z_crit_max - sqrt(w1_star) * z1_stage1) / sqrt(w2_star);
    val_min1 = min(val_11, val_12);
    alpha_c1_max = 1 - probnorm(val_min1);
    
    * calculate values within Phi function for H_{02};
    val_21 = (z_crit_max - sqrt(w1) * z2_stage1) / sqrt(w2);
    val_22 = (z_crit_max - sqrt(w1_star) * z2_stage1) / sqrt(w2_star);
    val_min2 = min(val_21, val_22);
    alpha_c2_max = 1 - probnorm(val_min2);
    
    t_crit1 = tinv((1 - alpha_c1_max), errdf_stage1);
    t_crit2 = tinv((1 - alpha_c2_max), errdf_stage1);
    
    keep z_: val_: alpha_: t_:;
run;

*--------------------------------------------------;
* estimated conditional target power;
*--------------------------------------------------;
data condpower;
    set achieved_power_stage1;
    planned_beta = 1 - &power_min.;
    beta1_hat = 1 - power; * power = power_stage1_alpha_max;
    cond_power = (beta1_hat - planned_beta) / beta1_hat;
    keep planned_beta beta1_hat power cond_power;
run;

*--------------------------------------------------;
* adaptive planning;
*--------------------------------------------------;
data adaptplan;
    set cl90pct1;
    adapt_planned_BE_ratio = &planned_BE_ratio.;
    adapt_planned_delta = log(adapt_planned_BE_ratio);
    delta_planned_ratio = log(&planned_BE_ratio.);
    if(sign(adapt_planned_delta) ^= sign(delta_hat_stage1)) then do;
        adapt_planned_BE_ratio = 1 / &planned_BE_ratio.;
    end;
    keep adapt_planned: delta_:;
run;


*-------------------------------------------------------------------;
* calculate sample size of second stage to achieve conditional power
*-------------------------------------------------------------------;
*copy variable to macro variable;
proc sql noprint;
    select adapt_planned_BE_ratio
    into  :adapt_planned_BE_ratio
    from adaptplan;
    
    select alpha_c1_max, alpha_c2_max 
    into  :alpha_c1_max,:alpha_c2_max
    from conderr;
    
    select var_hat_stage1
    into  :var_hat_stage1
    from iscv1;

    select cond_power
    into  :cond_power
    from condpower;
quit;

*--------------------------------------------------;
* power function for unequal significance levels;
*--------------------------------------------------;
%include ".\probmvt.sas"; *specify the folder;

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
        ABSEPS = .000001;
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

**********************************************************************;
* Final Analysis;
**********************************************************************;
*---------------------------------------;
* Analysis of second-stage data;
*---------------------------------------;
ods output Statistics=stat2 ConfLimits=cl2 EquivLimits=ev2 EquivTests=et2; 
proc ttest data=stage1plus2 tost(%sysfunc(log(&theta1.)), %sysfunc(log(&theta2.))) alpha=&alpha. order=internal;
    where STAGE = 2;
    var lnCmax1 lnCmax2 / crossover=(TRTA1 TRTA2) ;
run;

*---------------------------------------;
* check delta_hat_stage2 and 90pct CI;
*---------------------------------------;
data cl90pct2;
    set ev2;
    where Treatment = "Diff (1-2)" and Method="Pooled";
    delta_hat_stage2 = exp(-mean); * minus to convert diff(R - T) to diff(T - R);
    lower90_stage2 = exp(-upperCLMean);
    upper90_stage2 = exp(-lowerCLMean);
run;

*---------------------------------------;
*apply TOST at end of Stage 2;
*---------------------------------------;
data et2a;
    set et2;
    where Treatment = "Diff (1-2)" and Method="Pooled" and Test ^= "Overall";
    if Test = "上限" then Test = "Upper";
    if Test = "下限" then Test = "Lower";
run;

proc transpose data = et2a out = et2b prefix=P2nd_;
    by DF; * to keep in dataset;
    var probt;
    id Test;
run;

data tost2nd;
    set et2b;
    * convert p-values into z-statistics;
    z1_stage2 = probit(1 - P2nd_Lower);
    z2_stage2 = probit(1 - P2nd_Upper);
    errdf_stage2 = DF;
    keep p2nd_: z: errdf_stage2;
run;

*--------------------------------------------------;
* check Intra-subject CV of second stage;
*--------------------------------------------------;
data iscv2;
    set stat2(where=(Sequence ^= "Both" and Treatment = "Diff (1-2)")) end=last;
    retain _N _StdDev;
    
    if last then do;
        var_hat_stage2 = 0.5 * ((_N-1)*_StdDev**2 + (N-1)*StdDev**2) / (N+_N-2);
        cv_hat_stage2 = sqrt(exp(var_hat_stage2)-1);
        n_stage2 = _N + N;
        output;
     end;
    _N = N;
    _StdDev = StdDev;
run;

*--------------------------------------------------;
* Maximum Combination Test;
*--------------------------------------------------;
data maxcomb;
    merge tost1st tost2nd;
    w = &w.;
    w_star = &w_star.;
    z_crit_max = &z_crit_max.;
    
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
