/*###########################################################################
This is the example SAS code for Simulation of the Maximum Combination Test 
     proposed in the Research Article :

  Maurer, W, Jones, B, Chen, Y. 
  Controlling the type I error rate in two-stage sequential adaptive designs 
  when testing for average bioequivalence. 
  Statistics in Medicine. 2018; 37: 1587– 1607. 
  https://doi.org/10.1002/sim.7614

2021.10.30: Yusuke Morita, new for SAS User Group Meeting in Japan 2021;

############################################################################*/
/*****************************************************
configulation
*****************************************************/
%let w = 0.5;           * weight of maximum combination test ;
%let w_star = 0.25;     * weight of maximum combination test ;
%let alpha = 0.05;      * planned alpha (one-sided);
%let power = 0.80;      * planned power;
%let true_gmr = 0.95;   * planned geo-mean ratio;
%let theta1 = 0.8;      * BE lower limit (on ratio scale);
%let theta2 = 1.25;     * BE upper limit (on ratio scale);
%let futil_lower = 0.95;* lower limit of futility at stage 1 (on ratio scale);
%let futil_upper = 1.05;* upper limit of futility at stage 1 (on ratio scale);
%let n1 = 12;           * number of Subjects for stage 1;
%let n2_max = 100;      * maximum number of Subjects for stage 2;
%let sim_num = 300;     * repeat number of simulation per condition;

/*****************************************************
utility macros
*****************************************************/
%macro ODSOff(); /* call prior to BY-group processing */
    ods graphics off;
    ods exclude all; /* all open destinations */
    ods results off; /* no updates to tree view */
    options nonotes; /* optional, but sometimes useful */
%mend;

%macro ODSOn(); /* call after BY-group processing */
    ods graphics on;
    ods exclude none;
    ods results on;
    options notes;
%mend;

%macro calcAchievedPower;
    %let sigma2 = %sysfunc(log(1+&iscv.**2));* macro variable for the within-subject variance;
    %let std_derived = %sysevalf((&sigma2./2)**0.5); * macro variable for the derived common standard deviation;
    %let log_pt_8 = %sysfunc(log(&theta1.)); * macro variable for log(0.8);
    %let log_1_pt_25 = %sysfunc(log(&theta2.)); * macro variable for log(1.25);
    %let log_true_gmr = %sysfunc(log(&planned_BE_ratio.)); * macro variable for log(planned_BE_ratio);
    
    %ODSOff(); 
    ods output output = _achieved_power_stage1;
    ods select output;
    proc power;
         twosamplemeans test=equiv_diff alpha= &alpha_max.
         lower=&log_pt_8. upper=&log_1_pt_25. std=&std_derived.
         meandiff=&log_true_gmr.
         ntotal= &n_stage1.
         power = .;
    run;
    
    proc sql noprint;
        select power format=best12. into :power
        from _achieved_power_stage1;
    quit;
    %ODSOn(); 
    
%mend;

/*****************************************************
utility functions
*****************************************************/
proc fcmp outlib=WORK.FUNCDT.CARCALC;
    /* power by non-central t approximation (fast) */
    function power_nct(_alpha, _true_ratio, _lnstddev, _n , _df) ;
        _t1 = tinv(1 - _alpha, _df);
        _t2 = - _t1;
        _nc1 = sqrt(_n) * log(_true_ratio / 0.80) / sqrt(2) / _lnstddev;
        _nc2 = sqrt(_n) * log(_true_ratio / 1.25) / sqrt(2) / _lnstddev;
        _power = probt(_t2, _df, _nc2) - probt(_t1, _df, _nc1);
        return (_power) ;
    endsub ;

   /* exact power by proc power */
   function calcAchievedPower(alpha_max, n_stage1, iscv, planned_BE_ratio, theta1, theta2);
      rc=run_macro('calcAchievedPower', alpha_max, n_stage1, iscv, planned_BE_ratio, theta1, theta2, power);
      if rc = 0 then return(power);
      else return(.);
   endsub;
   
run ;

options cmplib=WORK.FUNCDT ; *--- specify stored library of custom functions ;



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


*****************************************************
Simulatoin data for stage 1
*****************************************************;
data stage1;
    call streaminit(20211118);
    stage = 1;
    do ratio = 0.80, 0.95;
        do iscv = 0.2 to 0.4 by 0.1;
            lnsigma = sqrt(log(1+iscv**2));            
            do N1 = &n1.;
                do trial = 1 to &sim_num.;
                    do subject = 1 to N1;
                        lnt = rand('normal', 1, sqrt(log(1+0.5**2))); /* 0.5 : InterSubject Variability , no effect on simulation */
                        lndif = rand('normal', log(ratio), lnsigma * sqrt(2));/* Var(ln(T/R)) = Var(ln(T)) + Var(ln(R)) = 2σ^2 */
                        lnr = lnt - lndif;
                        if mod(subject, 2) = 0 then do;
                            TRTSEQA = "TR";
                            TRTA1 = "T";
                            TRTA2 = "R";
                            lnAVAL1 = lnt;
                            lnAVAL2 = lnr;
                        end;
                        else do; 
                            TRTSEQA = "RT";
                            TRTA1 = "R";
                            TRTA2 = "T";
                            lnAVAL1 = lnr;
                            lnAVAL2 = lnt;
                        end;
                        output;
                    end;
                end;
            end;
        end;
    end;
run;

/*** check T/R ratio ********;
proc means data = stege1 n mean;
    by ratio iscv N;
    var  lnt lnr;
    output out = conf1 mean = meant meanr;
run;

data conf2;
    set conf1;
    t_r = exp(meant - meanr);
run;
*****************************/

**********************************************************************;
* Interim Analysis;
**********************************************************************;
%ODSOff(); 
ods output Statistics=stat1 ConfLimits=cl1 EquivLimits=ev1 EquivTests=et1; 
proc ttest data=stage1 tost(%sysfunc(log(&theta1.)), %sysfunc(log(&theta2.))) alpha=&alpha. order=internal;
    by ratio iscv n1 trial ;
    var lnAVAL1 lnAVAL2 / crossover=(TRTA1 TRTA2) ;
run;
%ODSOn(); 

*---------------------------------------;
* check futility;
*---------------------------------------;
data cl90pct1;
    set ev1;
    where Treatment = "Diff (1-2)" and Method="Pooled";
    delta_hat_stage1 = exp(-mean); * minus to convert diff(R - T) to diff(T - R);
    lower90_stage1 = exp(-upperCLMean);
    upper90_stage1 = exp(-lowerCLMean);
    
    if (upper90_stage1 < &futil_lower.) or (&futil_upper. < lower90_stage1) then futility1 = 1;
    else futility1 = 0;
    keep ratio iscv n1 trial delta_hat_stage1 lower90_stage1 upper90_stage1 futility1;
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
    by ratio iscv n1 trial DF; * to keep in dataset;
    var probt;
    id Test;
run;

data tost1st;
    merge et1b cl90pct1;
    by ratio iscv n1 trial;
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
    z_nABE_passed_stage1 = (z_nrejlow_stage1 & z_nrejupp_stage1 & futility1 = 0); *1:ABE accepted, 0:ABE rejected;
    
    keep ratio iscv n1 trial p1st_: z: DF futility1;
run;

*--------------------------------------------------;
* check achieved power of first stage1 for futility;
*--------------------------------------------------;
data acpower1;
    set stat1(where=(Sequence ^= "Both" and Treatment = "Diff (1-2)"));
    by ratio iscv n1 trial;
    retain _N _StdDev;
    if first.trial then do;
        _N = .;
        _StdDev = .;
    end;
    
    if last.trial then do;
        var_hat_stage1 = 0.5 * ((_N-1)*_StdDev**2 + (N-1)*StdDev**2) / (N+_N-2);
        cv_hat_stage1 = sqrt(exp(var_hat_stage1)-1);
        n_stage1 = _N + N;
        power1_nct = power_nct(&alpha_max., &true_gmr., sqrt(var_hat_stage1), n_stage1, n_stage1 - 2 );
        power1_exact = calcAchievedPower(&alpha_max., n_stage1, cv_hat_stage1, &true_gmr., &theta1., &theta2.);
        output;
     end;
    _N = N;
    _StdDev = StdDev;
    
    keep ratio iscv n1 trial var_hat_stage1 cv_hat_stage1 n_stage1 power1:;
run;

*********************************************************************;
* Sample Size Re-Estimation;
*********************************************************************;
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
    
    keep ratio iscv n1 trial z_: val_: alpha_: t_:;
run;

*--------------------------------------------------;
* estimated conditional target power;
*--------------------------------------------------;
data condpower;
    set acpower1;
    planned_beta = 1 - &power.;
    beta1_hat = 1 - power1_exact;*achieved power;
    cond_power = (beta1_hat - planned_beta) / beta1_hat;
    keep ratio iscv n1 trial planned_beta beta1_hat power: cond_power;
run;

*--------------------------------------------------;
* adaptive planning;
*--------------------------------------------------;
data adaptplan;
    set cl90pct1;
    adapt_planned_BE_ratio = &true_gmr.;
    adapt_planned_delta = log(adapt_planned_BE_ratio);
    delta_planned_ratio = log(&true_gmr.);
    if(sign(adapt_planned_delta) ^= sign(delta_hat_stage1)) then do;
        adapt_planned_BE_ratio = 1 / &true_gmr.;
    end;
    keep ratio iscv n1 trial adapt_planned: delta_:;
run;

*-------------------------------------------------------------------;
* calculate sample size of second stage to achieve conditional power
*-------------------------------------------------------------------;
*input for sample size re-estimation;
data inputssr;
    merge
        adaptplan
        conderr
        condpower
        acpower1
        cl90pct1
    ;
    by ratio iscv n1 trial;
    keep ratio iscv n1 trial adapt_planned_BE_ratio alpha_c1_max alpha_c2_max var_hat_stage1 cv_hat_stage1 cond_power futility1;
run;

*--------------------------------------------------;
* power function for unequal significance levels;
*--------------------------------------------------;
%include ".\probmvt.sas"; *specify the folder;

proc iml;
    load ; *load modules in probmvt.sas;

    varNames = {ratio iscv n1 trial adapt_planned_BE_ratio alpha_c1_max alpha_c2_max var_hat_stage1 cond_power};
    use inputssr(where=(futility1=0));
    ssr = {. . . . . . . .}; /** 1x8 numerical vector of results for a trial **/
    create outputssr from ssr[colname={"ratio" "iscv" "n1" "trial" "var_hat_stage1" "cond_power" "condpower" "n2"}];
 
    setin inputssr; /** make current for reading **/
    setout outputssr; /** make current for writing **/
 
    do data;
        read next var varNames;
        lnGmrPlanned = log(adapt_planned_BE_ratio);
        ln1p25 = log(&theta2.);
        variance1 = var_hat_stage1;
        
        COVAR = {1 -1, 
                -1  1};
        acheived = 0;
        
        do N = 2 to &n2_max. by 2 until(acheived);
            DIM = 2;
            DF = N - 2;
            if (DF = 0) then DF = 1; *when DF = 0 calculation fails;
            
            nc1 = (ln1p25 + lnGmrPlanned) /sqrt(variance1 * 2 / N);
            nc2 = (ln1p25 - lnGmrPlanned) /sqrt(variance1 * 2 / N); 
            DELTA = j(1, 2);
            DELTA[1] = nc1;
            DELTA[2] = nc2;
        
            INFIN = J(1, DIM, 1);
            tc1 = tinv(1 - alpha_c1_max, DF);
            tc2 = tinv(1 - alpha_c2_max, DF);
            LOWER = j(1, 2);
            LOWER[1] = tc1;
            LOWER[2] = tc2;     
            UPPER = j(1, 2, 999);*infinity:999 is ignored;
            
            MAXPTS = 2000*DIM*DIM*DIM;
            ABSEPS = .000001;
            RELEPS = 0;
      
            RUN MVN_DIST( DIM, DF, DELTA, LOWER, UPPER, INFIN, COVAR, MAXPTS, ABSEPS, RELEPS,  ERROR, condpower, NEVALS, INFORM );
    
            if (condpower >= cond_power) then acheived = 1;
            *print N condpower;
            if acheived = 1 then do;
                ssr[1]=ratio;
                ssr[2]=iscv;
                ssr[3]=n1;
                ssr[4]=trial;
                ssr[5]=var_hat_stage1;
                ssr[6]=cond_power;
                ssr[7]=condpower;
                ssr[8]=N;
                append from ssr;  
                goto breakout;
            end;
        end;
        *print N;
        breakout:
    end;
    close inputssr outputssr;
quit;

*--------------------------------------------------;
* integration of first stage results;
*--------------------------------------------------;
data interim;
    merge 
        tost1st
        acpower1
        conderr
        condpower
        adaptplan
        outputssr
    ;
    by ratio iscv n1 trial;
    w = &w.;
    w_star = &w_star.;
    futility2 = (futility1 = 0 and z_nABE_passed_stage1 = 0 and power1_exact >= &power.);
    futility3 = (futility1 = 0 and z_nABE_passed_stage1 = 0 and power1_exact <  &power. and n2 < 4);
    go_stage2 = (z_nABE_passed_stage1 = 0 and power1_exact <  &power. and n2 >= 4);
    
    *_check = sum(futility1,z_nABE_passed_stage1,futility2,futility3,go_stage2);
run;


**********************************************************************;
* Final Analysis;
**********************************************************************;
*****************************************************
Simulatoin data for stage 2
*****************************************************;
data stage2;
    call streaminit(20211118);
    set interim(where=(go_stage2 = 1) keep = ratio iscv n1 trial n2 go_stage2);
    stage = 2;
    lnsigma = sqrt(log(1+iscv**2));            
    do subject = 2001 to 2000 + N2;
        lnt = rand('normal', 1, sqrt(log(1+0.5**2))); /* 0.5 : InterSubject Variability , no effect on simulation */
        lndif = rand('normal', log(ratio), lnsigma * sqrt(2));/* Var(ln(T/R)) = Var(ln(T)) + Var(ln(R)) = 2σ^2 */
        lnr = lnt - lndif;
        if mod(subject, 2) = 0 then do;
            TRTSEQA = "TR";
            TRTA1 = "T";
            TRTA2 = "R";
            lnAVAL1 = lnt;
            lnAVAL2 = lnr;
        end;
        else do; 
            TRTSEQA = "RT";
            TRTA1 = "R";
            TRTA2 = "T";
            lnAVAL1 = lnr;
            lnAVAL2 = lnt;
        end;            
        output;
    end;
run;

*---------------------------------------;
* Analysis of second-stage data;
*---------------------------------------;
%ODSOff(); 
ods output Statistics=stat2 ConfLimits=cl2 EquivLimits=ev2 EquivTests=et2; 
proc ttest data=stage2 tost(%sysfunc(log(&theta1.)), %sysfunc(log(&theta2.))) alpha=&alpha. order=internal;
    by ratio iscv n1 trial ;
    var lnAVAL1 lnAVAL2 / crossover=(TRTA1 TRTA2) ;
run;
%ODSOn(); 

*---------------------------------------;
* apply TOST at end of Stage 2;
*---------------------------------------;
data et2a;
    set et2;
    where Treatment = "Diff (1-2)" and Method="Pooled" and Test ^= "Overall";
    if Test = "上限" then Test = "Upper";
    if Test = "下限" then Test = "Lower";
run;

proc transpose data = et2a out = et2b prefix=P2nd_;
    by ratio iscv n1 trial DF; * to keep in dataset;
    var probt;
    id Test;
run;

data tost2nd;
    set et2b;
    * convert p-values into z-statistics;
    z1_stage2 = probit(1 - P2nd_Lower);
    z2_stage2 = probit(1 - P2nd_Upper);
    errdf_stage2 = DF;
    keep ratio iscv n1 trial p2nd_: z: errdf_stage2;
run;


*--------------------------------------------------;
* Maximum Combination Test;
*--------------------------------------------------;
data final;
    merge tost1st tost2nd(in=in2nd);
    by ratio iscv n1 trial;
    if in2nd;
    
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
    nrejlow_stage2 = (z1_max >= z_crit_max);
    
    * reject null hypothesis at right side?;
    nrejupp_stage2 = (z2_max >= z_crit_max);
    
    * reject null hypothesis at both sides?,;
    * i.e., is ABE accepted?;
    z_passed_stage2 = (nrejlow_stage2 & nrejupp_stage2); *1:ABE accepted, 0:ABE rejected;

    keep ratio iscv n1 trial z_crit_max z1_w z1_w_star z1_max z2_w z2_w_star z2_max nrej: z_passed:;
run;

**********************************************************************;
* Simulation Summary;
**********************************************************************;
data summary;
    merge
        interim(keep= ratio iscv n1 trial w w_star z_nABE_passed_stage1 futility: go_stage2 n2)
        final(  keep= ratio iscv n1 trial z_passed_stage2)
    ;
    by ratio iscv n1 trial ;
    
    n = sum(n1, n2);

    if missing(z_passed_stage2) then z_passed_stage2 = 0;

    futility = (sum(futility1, futility2, futility3) > 0);
    z_passed_any = (sum(z_nABE_passed_stage1, z_passed_stage2) > 0);
run;

proc means data = summary noprint;
    by ratio iscv n1 w w_star;
    var futility1 z_nABE_passed_stage1 futility2 futility3 futility go_stage2 z_passed_stage2 z_passed_any;
    output out = result1 mean=/autoname;
run;

proc means data = summary noprint;
    by ratio iscv n1;
    var n;
    output out = result2 mean=n_mean p5=n5 p50=n50 p80 = n80 p95=n95;
run;

data result;
    merge
        result1
        result2
    ;
    by ratio iscv n1;
    format n_mean 8.1;
    drop _TYPE_;
    alpha = &alpha.;
    power = &power.;
    true_gmr = &true_gmr.;
    futil_lower = &futil_lower.;
    futil_upper = &futil_upper.;
run;

**********************************************************************;
* Report;
**********************************************************************;
option orientation=landscape;
ods pdf file = "/home/u43325994/Max-Comb-Test-Simulation.pdf";

proc report data = result split="\";
    columns (_freq_ w w_star alpha power futil_lower futil_upper n1 
                ratio iscv z_nABE_passed_stage1_mean futility_mean go_stage2_mean 
                z_passed_stage2_mean z_passed_any_mean 
                n_mean n5 n50 n80 n95);
    define _freq_ / order "repeat";
    define w/order "w";
    define w_star/order "w*";
    define alpha/order "alpha";
    define power/order "planned\power";
    define futil_lower/order "lower";
    define futil_upper/order "upper";
    define n1/order "N stage1";
    define ratio/order order=internal "Ratio";
    define iscv/order order=internal format=8.2 "True ISCV";
    define z_nABE_passed_stage1_mean /"BE stage1";
    define futility_mean /"Futility";
    define go_stage2_mean /"Proceed\stage2";
    define z_passed_stage2_mean /"BE stage2";
    define z_passed_any_mean /"BE any";
    define n_mean /"N Average";
    define n5 /'N 5%';
    define n50 /'N 50%';
    define n80 /'N 80%';
    define n95 /'N 95%';
run;
quit;

ods pdf close;
