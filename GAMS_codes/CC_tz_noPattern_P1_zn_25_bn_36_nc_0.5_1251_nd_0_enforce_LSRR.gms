$inlinecom   /* */
$offlisting
$offsymxref offsymlist
option optcr = 0;
option solver=gurobi;
option lp=gurobi;
option ResLim=3600;

* Saving datos computational time
Scalar
   start 'date + time at start of GAMS job'
   now   'current date + time'
   elapsed ' run + time'
   totalRunTime 'total run + time';

* Problem information
* bn borders nodes
* zn bus stations
* be(bn,bn) district borders
* file name

set           bn, be(bn,bn), zn, ze(zn,zn), p, t, sp(zn,zn,zn,zn), zebe(zn,zn,bn,bn), g(zn,zn);
parameter     rev(zn,zn,p,t),maxRevP(p),ptt(p,zn,zn,t), frac(p,zn,zn,t), trips(zn,zn);
singleton set cbd, pp(p);
alias         (zn,h,i,j,k);
alias         (bn,u,v,n,m);
set           iLj(i,j), uLv(u,v), zzpt(i,j,p,t);
scalar        r, xdim, ydim;

$GDXIN  "output_zn_25_bn_36_nc_0.75_2055.gdx"
$LOADDC bn be zn ze p t rev sp zebe xdim ydim g trips frac ptt
$GDXIN

iLj(i,j)      = (ord(i) < ord(j));
uLv(u,v)      = (ord(u) < ord(v) and be(u,v));
zzpt(iLj,p,t) = (rev(iLj,p,t));

* Result variables and sets
variables        handle,bestW,handleTotalrev;
parameter ntz(i,j),zone;
sets Ycopy(v,u),Wcopy(v,u),check,jL,list;
scalar maxRevenue, accumualtive;
scalar bestP;
parameter zoneborder, lambda;
* Problem variables
variables        F;
binary variables X,Y;
Positive Variable W;




set zz(zn,zn), speb(i,i,i,i,u,u);
zz(iLj) =  yes;
speb(i,j,iLj(h,k),u,v) = ((sp(i,j,h,k) or sp(i,j,k,h)) and zebe(h,k,u,v));

equations maxrev, tariffborders, tariffstairs, border, flow,coupling1,coupling2;

maxrev                 .. F - sum(zzpt(zz(i,j),pp,t),
                          (rev(i,j,pp,t))*X(zzpt)) =e= 0;
tariffborders(zz) .. sum(zzpt(zz,pp,t), ord(t)*X(zzpt))
                     - sum(speb(zz,h,k,u,v), (Y(u,v)+Y(v,u)))              =e= 1;
tariffstairs(zz)  .. sum(zzpt(zz,pp,t), X(zzpt))                  =e= 1;
border(uLv(u,v))  .. Y(v,u) + Y(u,v)                     =l= 1;
flow(u)           .. sum(be(u,v), W(be)) - sum(be(v,u), W(be))    =e= 0;
*degree(u)         .. sum(be(u,v), W(be)) + sum(be(v,u), W(be))    =l= 2;
coupling1(u,v) ..card(be)*Y(u,v)=g=W(u,v);
coupling2(u,v) .. Y(u,v)=l=W(u,v);

model tariffzones /all/;
*model tariffzones /maxrev, tariffborders, tariffstairs, border, flow, degree/;

Parameter status(*);



* ==========================================================
* Experiment A1: Store exact solutions for each price system
* ==========================================================
alias (p,ps,pt);


Parameter
    YA1(u,v,p)             "Optimal tariff-border decisions for price system p"
    WA1(u,v,p)             "Optimal border-flow values for price system p"
    XA1(i,j,p,t)           "Optimal charged-zone-count selection for OD pair (i,j) under price system p"

    maxRevA1(p)            "Maximum expected revenue under price system p"
    deltaRev(p)            "Percentage revenue loss relative to the highest revenue across price systems"

    modelStatA1(p)         "GAMS model status for price system p"
    solveStatA1(p)         "GAMS solver status for price system p"
    bestBoundA1(p)         "Best objective bound reported by the solver for price system p"
    runtimeA1(p)           "Solver runtime in seconds for price system p"
    gapA1(p)               "Relative optimality gap in percent for price system p"

    revTransfer(ps,pt)     "Expected revenue under target price pt using the source-optimal zoning from ps"
    transferLossAbs(ps,pt) "Absolute revenue loss from retaining the source-optimal zoning under target price pt"
    transferLossPct(ps,pt) "Transfer revenue loss as a percentage of optimal revenue under target price pt"

    zonesp(zn,p)           "Zone assignment of station zn under the zoning obtained for price system p"
    lambdap(p)             "Total number of tariff zones under price system p";


* ==============


maxRevenue = 0;
accumualtive = 0;
totalRunTime = 0;

pp(p) = no;
loop(p,
pp(p) = ord(p);
option MIP = gurobi;
$onecho > gurobi.opt
threads 0
$offecho
*Saving computatinal time
start=jnow;
* Problem
tariffzones.OptFile = 1;
solve tariffzones max F using mip;

* Store solver information for every price
modelStatA1(p) = tariffzones.modelstat;
solveStatA1(p) = tariffzones.solvestat;
bestBoundA1(p) = tariffzones.objest;
runtimeA1(p)   = tariffzones.resusd;
    
* Experiment A1: Store exact optimum only after certification
if (tariffzones.modelstat = 1
    and tariffzones.solvestat = 1,

    YA1(u,v,p)   = Y.l(u,v);
    WA1(u,v,p)   = W.l(u,v);
    XA1(i,j,p,t) = X.l(i,j,p,t);
    maxRevA1(p)  = tariffzones.objval;

);
gapA1(p) = abs(tariffzones.objest - tariffzones.objval)*100
           / max(1e-10, abs(tariffzones.objval));

handle.l(p)= (F.l);
if(F.l>maxRevenue,
Ycopy(v,u)=0;
Wcopy(u,v)=0;
bestP=ord(p);
maxRevenue=F.l;
Wcopy(u,v)=W.l(u,v);
Ycopy(v,u)=Y.l(u,v);
ntz(iLj) = 0;
ntz(iLj) = 1 + sum(speb(iLj,h,k,u,v), Y.l(u,v));
status('objval') = tariffzones.objval;status('objest') = tariffzones.objest;status('gapsolver') = ((tariffzones.objest-tariffzones.objval)/tariffzones.objest);
);

elapsed = (jnow - start)*24*3600;
totalRunTime=elapsed+totalRunTime;


** selecting arcs
zoneborder(i,j)$(g(i,j) and sum(be(n,m)$(zebe(i,j,n,m) and (Y.l(n,m) or Y.l(m,n))),1 ) >0) = 1;

lambda = 1; check(i) = no;
jL(i)=yes;
loop(i$jL(i), 
  jL(i)    = no; 
  list(j)  = no;
  check(j) = no;
  list(i)  = yes; 
  while((card(list)>card(check)),
   loop(k$(list(k) and not check(k) ),
    check(k) = yes; 
    loop(j$(g(k,j) and
            jL(j)   and not
            list(j) and 
            (zoneborder(k,j) + zoneborder(j,k) = 0)
           ),
    list(j) = yes; 
    jL(j)   = no; 
    );
   );
  );
  loop(k$list(k),
   zone(k) = lambda;
  ); 
  lambda = lambda + 1;
);
zoneborder(i,j) =0;

zonesp(zn, p)=zone(zn);
lambdap(p)=lambda-1;

zone(zn)=0;
lambda=0;
jL(i)=no;
list(j)  = no;
check(j) = no;
);


*

* ====================================================
* A1: Revenue portability across all price systems
* ====================================================

Parameter
    validOpt(p) 
;

validOpt(p) =(modelStatA1(p) = 1 and solveStatA1(p) = 1);

* Evaluate all ordered price-system transitions
revTransfer(ps,pt)$(validOpt(ps) and validOpt(pt)) =
    sum((iLj(i,j),t),rev(i,j,pt,t)*XA1(i,j,ps,t));

transferLossAbs(ps,pt)$(validOpt(ps) and validOpt(pt)) =maxRevA1(pt) - revTransfer(ps,pt);

transferLossPct(ps,pt)$( validOpt(ps) and validOpt(pt)and maxRevA1(pt) > 0) =100*transferLossAbs(ps,pt)/maxRevA1(pt);

    
* Compute revenue differences.
maxRevenue = smax(p$validOpt(p), maxRevA1(p));

deltaRev(p)$(validOpt(p) and maxRevenue > 0) =
    100*(maxRevenue-maxRevA1(p))/maxRevenue;
    

* ==========================================================
* Experiment A2: OD-level mechanism analysis
* ==========================================================

Parameter
    price(p,t)                 "Trip fare under price system p and charged-zone count t"
    demand(i,j,p,t)            "Expected PT demand for bidirectional OD pair (i,j) under price system p and zone count t"
    tripVolume(i,j)            "Total potential trips in both directions between i and j"
    probPT(i,j,p,t)            "Trip-weighted PT choice probability for bidirectional OD pair (i,j)"

    directEffect(ps,pt)         "Aggregate revenue effect of changing from source price ps to target price pt while retaining source zoning"
    totalEffect(ps,pt)          "Aggregate revenue effect of changing both the price system and optimal zoning"
    directEffectPct(ps,pt)      "Direct aggregate revenue effect as percentage of source optimal revenue"
    totalEffectPct(ps,pt)       "Total aggregate revenue effect as percentage of source optimal revenue"

    deltaCDirect(i,j,ps,pt)     "Direct OD fare effect from changing the price system while retaining source zoning"
    deltaCSpatial(i,j,ps,pt)    "Spatial OD fare effect from adapting zoning under target price system"
    deltaCTotal(i,j,ps,pt)      "Total OD fare effect from changing both price system and zoning"

    deltaRDirect(i,j,ps,pt)     "Direct OD revenue effect from changing the price system while retaining source zoning"
    deltaRSpatial(i,j,ps,pt)    "Spatial OD revenue effect from adapting zoning under target price system"
    deltaRTotal(i,j,ps,pt)      "Total OD revenue effect from changing both price system and zoning"

    deltaTSpatial(i,j,ps,pt)    "Change in charged-zone count for OD pair (i,j) caused by zoning adaptation"

    deltaDDirect(i,j,ps,pt)     "Direct OD demand effect from changing the price system while retaining source zoning"
    deltaDSpatial(i,j,ps,pt)    "Spatial OD demand effect from adapting zoning under target price system"
    deltaDTotal(i,j,ps,pt)      "Total OD demand effect from changing both price system and zoning"
    deltaPDirect(i,j,ps,pt)   "Direct OD PT choice-probability effect"
    deltaPSpatial(i,j,ps,pt)  "Spatial OD PT choice-probability effect"
    deltaPTotal(i,j,ps,pt)   "Total OD PT choice-probability effect";

* from instance generator!!
price(p,t)  = (1 + 0.2*(ord(p)-1)) * ord(t);    
    
******

directEffect(ps,pt)$(validOpt(ps) and validOpt(pt)) =revTransfer(ps,pt) - maxRevA1(ps);
    
directEffectPct(ps,pt)$(validOpt(ps) and validOpt(pt)and maxRevA1(ps) > 0) = 100 * (revTransfer(ps,pt) - maxRevA1(ps)) / maxRevA1(ps);
        
totalEffect(ps,pt)$(validOpt(ps) and validOpt(pt)) = maxRevA1(pt) - maxRevA1(ps);

totalEffectPct(ps,pt)$( validOpt(ps) and validOpt(pt) and maxRevA1(ps) > 0) =100 * (maxRevA1(pt) - maxRevA1(ps))/ maxRevA1(ps);
    
* Fare

* Direct effect: change price, keep source zoning
deltaCDirect(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) =
  sum(t,price(pt,t)*XA1(i,j,ps,t))-sum(t,price(ps,t)*XA1(i,j,ps,t));

* Spatial effect: change zoning, keep target price
deltaCSpatial(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) =
  sum(t,price(pt,t)*XA1(i,j,pt,t))-sum(t,price(pt,t)*XA1(i,j,ps,t));

* Total effect: change both price and zoning -> optimalSol(fare(ps))- optimalSol(fare(pt)
deltaCTotal(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) =
    sum(t, price(pt,t)*XA1(i,j,pt,t))- sum(t, price(ps,t)*XA1(i,j,ps,t));
    
* REVENUE
* Direct effect: change price, keep source zoning
deltaRDirect(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) =
  sum(t,rev(i,j,pt,t)*XA1(i,j,ps,t))-sum(t,rev(i,j,ps,t)*XA1(i,j,ps,t));

* Spatial effect: change zoning, keep target price
deltaRSpatial(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) =
  sum(t,rev(i,j,pt,t)*XA1(i,j,pt,t))-sum(t,rev(i,j,pt,t)*XA1(i,j,ps,t));

* Total effect: change both price and zoning -> optimalSol(Rev(ps))- optimalSol(Revenue(pt))
deltaRTotal(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) =
    sum(t, rev(i,j,pt,t)*XA1(i,j,pt,t))- sum(t, rev(i,j,ps,t)*XA1(i,j,ps,t));

* zone counts

deltaTSpatial(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) =
sum(t, ord(t) *(XA1(i,j,pt,t) - XA1(i,j,ps,t)));


* Demand

demand(i,j,p,t)$iLj(i,j) =ptt(p,i,j,t)+ ptt(p,j,i,t);


* Direct demand effect
deltaDDirect(i,j,ps,pt)$( iLj(i,j) and validOpt(ps) and validOpt(pt)) =
    sum(t,(demand(i,j,pt,t)-demand(i,j,ps,t))*XA1(i,j,ps,t));

* Spatial demand effect
deltaDSpatial(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) = sum(t,demand(i,j,pt,t)*(XA1(i,j,pt,t)-XA1(i,j,ps,t)));
 
* Total demand effect optimalSol(Dem(ps))- optimalSol(Dem(pt))
deltaDTotal(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) = sum(t,demand(i,j,pt,t)*XA1(i,j,pt,t))-sum(t,demand(i,j,ps,t)*XA1(i,j,ps,t));

tripVolume(i,j)$iLj(i,j) = trips(i,j) + trips(j,i);

probPT(i,j,p,t)$(iLj(i,j) and tripVolume(i,j) > 0) =demand(i,j,p,t)/tripVolume(i,j);


* PT Share


* Direct effect: change price, keep source zoning
deltaPDirect(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) = sum(t,(probPT(i,j,pt,t) - probPT(i,j,ps,t))* XA1(i,j,ps,t));

* Spatial effect: change zoning, keep target price
deltaPSpatial(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) =sum(t,probPT(i,j,pt,t)* (XA1(i,j,pt,t) - XA1(i,j,ps,t)));

* Total effect: change price and zoning optimalSol(PT(ps))- optimalSol(PT(pt))
deltaPTotal(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) =sum(t,probPT(i,j,pt,t) * XA1(i,j,pt,t))- sum(t,probPT(i,j,ps,t) * XA1(i,j,ps,t));

* ================
* up to here A2
* ================

*=======
* * Mechanism of classification per OD Pair
*=======


Set cl "Categories of OD-level zoning effects" /
    Aplus   "More charged zones with revenue gain"
    Aminus  "More charged zones with revenue loss"
    Bplus   "Fewer charged zones with revenue gain"
    Bminus  "Fewer charged zones with revenue loss"
    C       "No change in charged-zone count"
    D       "Changed zone count with negligible revenue effect"
/

classification(i,j,ps,pt,cl) "Classification of each OD pair and price transition";

Scalar epsR /1e-7/;

classification(i,j,ps,pt,'Aplus')$(
    iLj(i,j)
    and validOpt(ps) and validOpt(pt)
    and ord(ps) <> ord(pt)
    and deltaTSpatial(i,j,ps,pt)>0
    and deltaRSpatial(i,j,ps,pt)>epsR
)=yes;

classification(i,j,ps,pt,'Aminus')$(
    iLj(i,j)
    and validOpt(ps) and validOpt(pt)
    and ord(ps) <> ord(pt)
    and deltaTSpatial(i,j,ps,pt)>0
    and deltaRSpatial(i,j,ps,pt)<-epsR
)=yes;

classification(i,j,ps,pt,'Bplus')$(
    iLj(i,j)
    and validOpt(ps) and validOpt(pt)
    and ord(ps) <> ord(pt)
    and deltaTSpatial(i,j,ps,pt)<0
    and deltaRSpatial(i,j,ps,pt)>epsR
)=yes;

classification(i,j,ps,pt,'Bminus')$(
    iLj(i,j)
    and validOpt(ps) and validOpt(pt)
    and ord(ps) <> ord(pt)
    and deltaTSpatial(i,j,ps,pt)<0
    and deltaRSpatial(i,j,ps,pt)<-epsR
)=yes;

classification(i,j,ps,pt,'C')$(
    iLj(i,j)
    and validOpt(ps) and validOpt(pt)
    and ord(ps) <> ord(pt)
    and deltaTSpatial(i,j,ps,pt)=0
)=yes;

classification(i,j,ps,pt,'D')$(
    iLj(i,j)
    and validOpt(ps) and validOpt(pt)
    and ord(ps) <> ord(pt)
    and abs(deltaTSpatial(i,j,ps,pt))>0
    and abs(deltaRSpatial(i,j,ps,pt))<=epsR
)=yes;
    


Parameter
    revIncreaseZones(ps,pt) "Net spatial revenue contribution from OD pairs with more charged zones"

    revDecreaseZones(ps,pt) "Net spatial revenue contribution from OD pairs with fewer charged zones"

    demandIncreaseZones(ps,pt) "Net spatial demand effect from OD pairs with more charged zones"

    demandDecreaseZones(ps,pt) "Net spatial demand effect from OD pairs with fewer charged zones"
    demandFrozen(i,j,ps,pt)   "PT demand under target price with source zoning"
    demandOptim(i,j,ps,pt)    "PT demand under target price with target zoning"
    fareFrozen(i,j,ps,pt)     "Fare under target price with source zoning"
    fareOptim(i,j,ps,pt)      "Fare under target price with target zoning"

    elPartA(i,j,ps,pt)        "Proportional spatial demand change"
    elPartB(i,j,ps,pt)        "Proportional spatial fare change"
    elasticityPT(i,j,ps,pt)   "Spatial arc elasticity of expected PT demand";


revIncreaseZones(ps,pt) =sum((i,j)$(iLj(i,j) and deltaTSpatial(i,j,ps,pt) > 0),deltaRSpatial(i,j,ps,pt));

revDecreaseZones(ps,pt) =sum((i,j)$(iLj(i,j)and deltaTSpatial(i,j,ps,pt) < 0),deltaRSpatial(i,j,ps,pt));

demandIncreaseZones(ps,pt) =sum((i,j)$(iLj(i,j) and deltaTSpatial(i,j,ps,pt) > 0),deltaDSpatial(i,j,ps,pt));

demandDecreaseZones(ps,pt) =sum((i,j)$(iLj(i,j) and deltaTSpatial(i,j,ps,pt) < 0), deltaDSpatial(i,j,ps,pt));


demandFrozen(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) = sum(t,demand(i,j,pt,t)*XA1(i,j,ps,t));

demandOptim(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) = sum(t,demand(i,j,pt,t)*XA1(i,j,pt,t));

fareFrozen(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) = sum(t,price(pt,t)*XA1(i,j,ps,t));

fareOptim(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt)) = sum(t,price(pt,t)*XA1(i,j,pt,t));
    

elPartA(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt) and demandFrozen(i,j,ps,pt)+demandOptim(i,j,ps,pt)>0) = 2*deltaDSpatial(i,j,ps,pt)/(demandFrozen(i,j,ps,pt)+demandOptim(i,j,ps,pt));

elPartB(i,j,ps,pt)$(iLj(i,j) and validOpt(ps) and validOpt(pt) and fareFrozen(i,j,ps,pt)+fareOptim(i,j,ps,pt)>0) = 2*deltaCSpatial(i,j,ps,pt) /(fareFrozen(i,j,ps,pt)+fareOptim(i,j,ps,pt));

elasticityPT(i,j,ps,pt)$( iLj(i,j) and validOpt(ps) and validOpt(pt) and abs(elPartB(i,j,ps,pt))>1e-10) = elPartA(i,j,ps,pt)/elPartB(i,j,ps,pt);

Parameter
    nOD(ps,pt,cl)       "Number of OD pairs in each category"
    revByClass(ps,pt,cl) "Net spatial revenue contribution by category"
    demByClass(ps,pt,cl) "Net spatial demand change by category";

nOD(ps,pt,cl) = sum((i,j)$classification(i,j,ps,pt,cl),1);

revByClass(ps,pt,cl) = sum((i,j)$classification(i,j,ps,pt,cl), deltaRSpatial(i,j,ps,pt));

demByClass(ps,pt,cl) = sum((i,j)$classification(i,j,ps,pt,cl), deltaDSpatial(i,j,ps,pt));

* ==========================================================
* A2: Consistency checks
* ==========================================================

Parameter
    checkClassification(ps,pt) "Difference between classified and total OD pairs"

    checkRevenue(ps,pt)  "Difference between class revenue contributions and A1 transfer loss"

    checkDemand(ps,pt) "Difference between class demand contributions and total spatial demand effect";

checkClassification(ps,pt)$(validOpt(ps) and validOpt(pt) and ord(ps) <> ord(pt)) = sum(cl,nOD(ps,pt,cl)) - sum((i,j)$iLj(i,j),1);

checkRevenue(ps,pt)$( validOpt(ps) and validOpt(pt) and ord(ps) <> ord(pt) ) = sum(cl,revByClass(ps,pt,cl)) - transferLossAbs(ps,pt);

checkDemand(ps,pt)$(validOpt(ps) and validOpt(pt) and ord(ps) <> ord(pt)) = sum(cl,demByClass(ps,pt,cl))- sum((i,j)$iLj(i,j),deltaDSpatial(i,j,ps,pt));


* ==========================================================
* A2: Additional mechanism summary statistics
* ==========================================================

Parameter
    shareOD(ps,pt,cl) "Percentage of OD pairs in each class"

    volumeByClass(ps,pt,cl) "Potential trip volume represented by each class"

    shareVolume(ps,pt,cl) "Percentage of potential trip volume in each class"

    revPositiveByClass(ps,pt,cl) "Positive spatial revenue contributions"

    revNegativeByClass(ps,pt,cl) "Negative spatial revenue contributions";

alias(cl,cll);
* Share of OD pairs
shareOD(ps,pt,cl)$(validOpt(ps) and validOpt(pt)and ord(ps) <> ord(pt) and sum(cll,nOD(ps,pt,cll)) > 0) = 100*nOD(ps,pt,cl) /sum(cll,nOD(ps,pt,cll));

* Potential trip volume by category
volumeByClass(ps,pt,cl) = sum((i,j)$classification(i,j,ps,pt,cl), tripVolume(i,j));

* Percentage of potential trips
shareVolume(ps,pt,cl)$( validOpt(ps) and validOpt(pt) and ord(ps) <> ord(pt) and sum(cll,volumeByClass(ps,pt,cll)) > 0 ) = 100*volumeByClass(ps,pt,cl) /sum(cll,volumeByClass(ps,pt,cll));

* Positive spatial revenue effects
revPositiveByClass(ps,pt,cl) = sum((i,j)$( classification(i,j,ps,pt,cl) and deltaRSpatial(i,j,ps,pt) > 0 ), deltaRSpatial(i,j,ps,pt));

* Negative spatial revenue effects
revNegativeByClass(ps,pt,cl) = sum((i,j)$(classification(i,j,ps,pt,cl) and deltaRSpatial(i,j,ps,pt) < 0), deltaRSpatial(i,j,ps,pt));











** selecting arcs
zoneborder(i,j)$(g(i,j) and sum(be(n,m)$(zebe(i,j,n,m) and (Ycopy(n,m) or Ycopy(m,n))),1 ) >0) = 1;
*zoneborder(i,j)$(g(i,j) and (sum(be(v,u)$zebe(i,j,u,v),Wcopy(u,v)) or sum(be(v,u)$zebe(i,j,v,u),Wcopy(v,u)) ))  = 1;



lambda = 1; check(i) = no;
jL(i)=yes;
loop(i$jL(i), 
  jL(i)    = no; 
  list(j)  = no;
  check(j) = no;
  list(i)  = yes; 
  while((card(list)>card(check)),
   loop(k$(list(k) and not check(k) ),
    check(k) = yes; 
    loop(j$(g(k,j) and
            jL(j)   and not
            list(j) and 
            (zoneborder(k,j) + zoneborder(j,k) = 0)
           ),
    list(j) = yes; 
    jL(j)   = no; 
    );
   );
  );
  loop(k$list(k),
   zone(k) = lambda;
  ); 
  lambda = lambda + 1;
);
zoneborder(i,j) =0;



********

scalar scale,xc,yc;
scale = 12.25/(xdim-1);

scalar noz, r,gg,b;
parameter first(zn), second(zn),third(zn);

noz = smax(zn, zone(zn) );

while (noz > 0,
 r = uniform(0,1);
 gg = uniform(0,1);
 b = uniform(0,1);
 loop (zn$(zone(zn) = noz),
 first(zn)=r;
 second(zn)=gg;
 third(zn)=b;
 );
 noz = noz - 1;
);


file tzres /tzresult_noPatten_nd_0_zn_25_bn_36_nc_0.75_2055_p1_A1Exp.tex/; 
put  tzres;

put '\documentclass{article}' /;
put '\usepackage{tikz}' /;
put '\usetikzlibrary{arrows.meta}' /;
put '\begin{document}' /;
put '\begin{tikzpicture}[scale = ' scale:7:5 ']' /;
yc = 0;
loop(bn,
 xc = mod(ord(bn),xdim +1);
 if (xc = 0, xc = xdim +1);
 if (xc = 1, yc = yc + 1);
 put '\coordinate (' bn.tl:5:0 ') at ('xc:5:0',' yc:5:0 ');' /;
);
yc = 0;
loop(zn,
 xc = mod(ord(zn),xdim);
 if (xc = 0, xc = xdim);
 if (xc = 1, yc = yc + 1);
** this is new here
put '\definecolor{currentcolor}{rgb}{'(first(zn))', '(second(zn))', '(third(zn))'};' /;  
 put '\node (' zn.tl:4:0 ') at ('(xc+0.5):4:1',' (yc+0.5):4:1 ') [fill=currentcolor, inner sep=1.5pt, draw, circle, text=blue, minimum height = 1cm]{{\tiny ' zn.tl:4:0'}};'  /;
);
loop(be(u,v)$(W.up(be)>0),
  put '\draw[arrows = -{Stealth[inset=0pt]}] (' u.tl:5:0')--('v.tl:5:0');'/;
);
loop(g(i,j),
  put '\draw[dashed,line width = 0.1cm,blue] (' i.tl:5:0')--('j.tl:5:0');'/;
);
loop(be(u,v)$(Wcopy(u,v) or Wcopy(v,u)),
put '\draw[inner sep=0.5pt,line width = 0.1cm,red] (' u.tl:5:0')--('v.tl:5:0');'/;

*  put '\draw[-,red] (' u.tl:5:0')--('v.tl:5:0 ');'/;
);
loop(bn,
  put '  \node at (' bn.tl:5:0 ') [inner sep=1.5pt, text=black, fill=white,draw]{{\tiny ' bn.tl:5:0'}};'/;
);
put '\end{tikzpicture}' /;
put '\end{document}' /;
putclose;

execute_unload 'tmp_noPattern_p1_nd_0_zn_25_bn_36_nc_0.75_2055_A1Exp.gdx';  

execute_unload 'exactOptima_noPattern_p1_nd_0_zn_25_bn_36_nc_0.75_2055_A1Exp.gdx',  
    bestP,
    maxRevA1,
    bestBoundA1,
    gapA1,
    runtimeA1,
    modelStatA1,
    solveStatA1,
    validOpt,
    lambdap;
    
execute_unload 'transitionMatrixExact_noPattern_p1_nd_0_zn_25_bn_36_nc_0.75_2055_A1Exp.gdx',  
    bestP,
    validOpt,
    maxRevA1,
    revTransfer,
    transferLossAbs,
    transferLossPct,
    directEffect,
    directEffectPct,
    totalEffect,
    totalEffectPct;
    
    
execute_unload 'ODlevelExact_noPattern_p1_nd_0_zn_25_bn_36_nc_0.75_2055_A2Exp.gdx',  
    bestP,
    zn,
    validOpt,
    iLj,
    XA1,
    price,
    demand,
    tripVolume,
    probPT,
    deltaTSpatial,
    deltaCDirect,
    deltaCSpatial,
    deltaCTotal,
    deltaDDirect,
    deltaDSpatial,
    deltaDTotal,
    deltaPDirect,
    deltaPSpatial,
    deltaPTotal,
    deltaRDirect,
    deltaRSpatial,
    deltaRTotal,
    demandFrozen,
    demandOptim,
    fareFrozen,
    fareOptim,
    elPartA,
    elPartB,
    elasticityPT,
    classification;
    
execute_unload 'MechanismSummary_noPattern_p1_nd_0_zn_25_bn_36_nc_0.75_2055_A2Exp.gdx',  
    price,
    cl,
    validOpt,
    nOD,
    shareOD,
    volumeByClass,
    shareVolume,
    revByClass,
    revPositiveByClass,
    revNegativeByClass,
    demByClass,
    revIncreaseZones,
    revDecreaseZones,
    demandIncreaseZones,
    demandDecreaseZones,
    transferLossAbs,
    checkClassification,
    checkRevenue,
    checkDemand;

$exit
