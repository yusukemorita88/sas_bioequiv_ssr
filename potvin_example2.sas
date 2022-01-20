/*###########################################################################
This is the example SAS code of the Research Article :

 Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ, Smith RA.
 Sequential design approaches for bioequivalence studies with crossover designs. 
 Pharm Stat. 2008 Oct-Dec;7(4):245-62. 
 doi: 10.1002/pst.294. PMID: 17710740.

2021.10.30 Yusuke Morita, new for SAS User Group Meeting in Japan 2021;

############################################################################*/
data example2;
input usubjid $ stage seq $ t r;
cards;
101 1 TR 4.8 4.6
102 1 TR 1.4 1.3
103 1 TR 3.9 4.4 
104 1 TR 1.5 2.4
105 1 TR 3.8 2.9
106 1 TR 3.3 2.5
107 1 RT 1.5 1.7
108 1 RT 1.9 1.3
109 1 RT 2.4 2.1
110 1 RT 13.9 11.4
111 1 RT 2.3 1.5
112 1 RT 2.4 2.6
213 2 TR 4.0 4.2
214 2 TR 0.7 0.5
215 2 TR 10.3 11.9
216 2 TR 4.7 4.8
217 2 RT 5.5 5.9
218 2 RT 5.4 3.8
219 2 RT 2.1 4.3
220 2 RT 2.4 3.6
;

run;

proc transpose data = example2 out = test2;
    by usubjid stage seq;
    var t r;
run;

data test3;
    set test2;
    by usubjid ;
    period = (substr(seq,1,1) ^= upcase(_name_));
    lnaval = log(col1);
    rename
     col1 = aval
     _name_ = trtp
    ;
run;

**********************************************************************;
* Method B;
**********************************************************************;

**** STAGE 1 *******************;
ods output LSMeanDiffCL = lsm1 OverallANOVA = anova1;
proc glm data = test3;
    where stage = 1;
    class trtp(ref='r') period usubjid seq;
    model lnaval = trtp period seq usubjid ;
    lsmeans trtp / cl diff=control('r') alpha=0.10 ;
run;
quit;

data stage1;
    merge
        anova1(where=(Source="Error"))
        lsm1
    ;
    estimate = exp(difference);
    lnlower = difference - tinv(1 - 0.0294, DF) * sqrt(MS * 2 / (DF+2));
    lowerp = exp(lnlower);
    lnupper = difference + tinv(1 - 0.0294, DF) * sqrt(MS * 2 / (DF+2));
    upperp = exp(lnupper);

run;

*Power;
data power;
    set anova1(where=(Source="Error"));
    n = DF + 2;
    alpha = 0.0294;
    *alpha = 0.05;
    intraSubjectV = sqrt(MS);
    uRatio = 0.95;
    n2 = n - 2;
    t1 = tinv(1 - alpha, n - 2);
    t2 = - t1;
    nc1 = sqrt(n)*log(uRatio / 0.80) / sqrt(2) / intraSubjectV;
    nc2 = sqrt(n)*log(uRatio / 1.25) / sqrt(2) / intraSubjectV;
    
    /*central-t approximation*/
    power = probt(t2 - nc2, df) - probt(t1 - nc1, df);

    /*non-central-t approximation*/
    *power = probt(t2, df, nc2) - probt(t1, df, nc1);
run;


data ssr;
    set anova1(where=(Source="Error"));
    n = (DF + 2) + 8; *stage1 + 8 subjects;
    alpha = 0.0294;
    intraSubjectV = sqrt(MS);
    uRatio = 0.95;
    df = n - 3;

    t1 = tinv(1 - alpha, df);
    t2 = - t1;
    nc1 = sqrt(n)*log(uRatio / 0.80) / sqrt(2) / intraSubjectV;
    nc2 = sqrt(n)*log(uRatio / 1.25) / sqrt(2) / intraSubjectV;

    /*central-t approximation*/
    power = probt(t2 - nc2, df) - probt(t1 - nc1, df);

    /*non-central-t approximation*/
    *power = probt(t2, df, nc2) - probt(t1, df, nc1);
run;


**** STAGE 1+2 *******************;
ods output LSMeanDiffCL = lsm2 OverallANOVA = anova2;
proc glm data = test3;
    where stage in (1, 2);
    class trtp(ref='r') period usubjid seq stage;
    model lnaval = trtp seq stage period(stage) usubjid/ ;
    lsmeans trtp /cl diff=control('r') alpha=0.10 ;
run;
quit;

data stage2;
    merge
        anova2(where=(Source="Error"))
        lsm2
    ;
    estimate = exp(difference);
    lnlower = Difference - tinv(1 - 0.0294, DF) * sqrt(MS * 2 / (DF+3));
    lowerp = exp(lnlower);
    lnupper = Difference + tinv(1 - 0.0294, DF) * sqrt(MS * 2 / (DF+3));
    upperp = exp(lnupper);

run;

**********************************************************************;
* Method C;
**********************************************************************;

**** STAGE 1 *******************;
*Power;
data powerC;
    set anova1(where=(Source="Error"));
    n = DF + 2;
    alpha = 0.05;
    intraSubjectV = sqrt(MS);
    uRatio = 0.95;
    n2 = n - 2;
    t1 = tinv(1 - alpha, n - 2);
    t2 = - t1;
    nc1 = sqrt(n)*log(uRatio / 0.80) / sqrt(2) / intraSubjectV;
    nc2 = sqrt(n)*log(uRatio / 1.25) / sqrt(2) / intraSubjectV;

    /*central-t approximation*/
    df = n - 2;
    power = probt(t2 - nc2, df) - probt(t1 - nc1, df);

    /*non-central-t approximation*/
    *power = probt(t2, df, nc2) - probt(t1, df, nc1);

run;
