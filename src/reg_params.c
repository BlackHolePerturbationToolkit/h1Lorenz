#include <math.h>
#include <gsl/gsl_sf_hyperg.h>
#include "math_funcs.h"



#define Power(x, y)    				 				(pow((double)(x), (double)(y)))
#define Sqrt(x)        				 				(sqrt((double)(x)))
#define EllipticK(x)									(elliptic_integral_1st((double)(sqrt(x))))
#define EllipticE(x)									(elliptic_integral_2nd((double)(sqrt(x))))
#define Pi														M_PI
#define Hypergeometric2F1(a,b,c,x)		(hypergeom2F1(a,b,c,x))

double M = 1.0;

double hypergeom2F1(double a, double b, double c, double x)
{
	return	gsl_sf_hyperg_2F1(a,b,c,x);
}

double Anna_and_Barrys_D_r(double l, double p, double e, double r1, double L)
{
	double e1 = Sqrt(((-2 - 2*e + p)*(-2 + 2*e + p))/(p*(-3 - Power(e,2) + p)));

	return -(((188*Power(L,10)*(2 - r1) + 46*Power(L,8)*Power(r1,2)*(34 + (-17 + 6*Power(e1,2))*r1) + 
        Power(r1,11)*(30 + 244*Power(e1,2) + (-15 - 14*Power(e1,2) + 30*Power(e1,4))*r1) + 
        Power(L,6)*Power(r1,4)*(2444 + 4*(-299 + 280*Power(e1,2))*r1 - (13 + 38*Power(e1,2))*Power(r1,2)) + 
        Power(L,4)*Power(r1,6)*(1700 + 8*(-96 + 211*Power(e1,2))*r1 + (-41 - 106*Power(e1,2) + 22*Power(e1,4))*Power(r1,2)) + 
        Power(L,2)*Power(r1,8)*(444 + 136*(-1 + 8*Power(e1,2))*r1 + (-43 - 82*Power(e1,2) + 68*Power(e1,4))*Power(r1,2)))*
      EllipticE(Power(L,2)/(Power(L,2) + Power(r1,2))))/
    (8.*(-0.5 + l)*(1.5 + l)*Pi*Power(r1,7)*Power(Power(L,2) + Power(r1,2),3.5)) + 
   ((-94*Power(L,8)*(2 - r1) - Power(r1,9)*(2*(15 + 92*Power(e1,2)) + (-15 + Power(e1,2) + 15*Power(e1,4))*r1) + 
        2*Power(L,6)*Power(r1,2)*(-284 + (128 - 105*Power(e1,2))*r1 + (7 + 6*Power(e1,2))*Power(r1,2)) + 
        Power(L,2)*Power(r1,6)*(-192 + (8 - 594*Power(e1,2))*r1 + (44 + 18*Power(e1,2) - 35*Power(e1,4))*Power(r1,2)) + 
        Power(L,4)*Power(r1,4)*(-572 + 20*(10 - 31*Power(e1,2))*r1 + (43 + 31*Power(e1,2) - 12*Power(e1,4))*Power(r1,2)))*
      EllipticK(Power(L,2)/(Power(L,2) + Power(r1,2))))/
    (8.*(-0.5 + l)*(1.5 + l)*Pi*Power(r1,5)*Power(Power(L,2) + Power(r1,2),3.5)));

}

double Anna_and_Barrys_D_t(double l, double p, double e, double r1, double ur1, double L)
{
	double e1 = Sqrt(((-2 - 2*e + p)*(-2 + 2*e + p))/(p*(-3 - Power(e,2) + p)));

	return -(((e1*(276*Power(L,8) + 1044*Power(L,6)*Power(r1,2) + 
          Power(L,4)*Power(r1,4)*(1550 + (-37 + 22*Power(e1,2))*r1) + 
          Power(r1,8)*(274 + (-29 + 30*Power(e1,2))*r1) + 
          2*Power(L,2)*Power(r1,6)*(528 + (-33 + 34*Power(e1,2))*r1))*ur1*
        EllipticE(Power(L,2)/(Power(L,2) + Power(r1,2))))/
      (8.*Pi*(-2 + r1)*Power(r1,3)*Power(Power(L,2) + Power(r1,2),3.5)) - 
     (e1*(186*Power(L,6) + 6*Power(L,4)*Power(r1,2)*(99 + (-3 + 2*Power(e1,2))*r1) + 
          Power(r1,6)*(214 + (-14 + 15*Power(e1,2))*r1) + Power(L,2)*Power(r1,4)*(622 + (-32 + 35*Power(e1,2))*r1)
          )*ur1*EllipticK(Power(L,2)/(Power(L,2) + Power(r1,2))))/
      (8.*Pi*(-2 + r1)*r1*Power(Power(L,2) + Power(r1,2),3.5)))/((-0.5 + l)*(1.5 + l)));
}

double Anna_and_Barrys_D_H(double l, double r0, double E0, double L0)
{
	return (-((Power(L0,4)*(-36 + 8*r0) + Power(L0,2)*Power(r0,2)*(-74 + 17*r0) + Power(r0,4)*(-38 + (9 - 2*Power(E0,2))*r0))*
         EllipticE(Power(L0,2)/(Power(L0,2) + Power(r0,2))))/(4.*Pi*Power(r0,3)*Power(Power(L0,2) + Power(r0,2),1.5))\
      - ((8*Power(L0,4)*(4 - r0) + 2*Power(L0,2)*Power(r0,2)*(49 + 8*(-2 + Power(E0,2))*r0) + 
          Power(r0,4)*(66 + (-24 + 17*Power(E0,2))*r0))*EllipticK(Power(L0,2)/(Power(L0,2) + Power(r0,2))))/
      (4.*Pi*Power(r0,3)*Power(Power(L0,2) + Power(r0,2),1.5)))/((-0.5 + l)*(1.5 + l));
}

double Anna_and_Barrys_F_H(double l, double r0, double E0, double L0)
{
	return (3*(1536*Power(L0,14)*Power(M,2) + 320*Power(L0,12)*M*(53*M - 5*r0)*Power(r0,2) + 8*Power(L0,10)*M*Power(r0,4)*(8737*M + 4*(-395 + 112*Power(E0,2))*r0) + 
        4*Power(L0,8)*Power(r0,6)*(35486*Power(M,2) + 7*(-1235 + 582*Power(E0,2))*M*r0 + 30*Power(r0,2)) + 
        Power(L0,6)*Power(r0,8)*(153296*Power(M,2) + 2*(-21495 + 14114*Power(E0,2))*M*r0 + 15*(17 + 8*Power(E0,2))*Power(r0,2)) + 
        Power(L0,4)*Power(r0,10)*(85100*Power(M,2) + 2*(-11975 + 10849*Power(E0,2))*M*r0 + 15*(4 + 57*Power(E0,2) - 32*Power(E0,4))*Power(r0,2)) + 
        Power(r0,14)*(340*Power(M,2) + 34*(25 + Power(E0,2))*M*r0 - 15*(6 - 33*Power(E0,2) + 28*Power(E0,4))*Power(r0,2)) + 
        Power(L0,2)*Power(r0,12)*(19616*Power(M,2) + 2*(-1825 + 3108*Power(E0,2))*M*r0 - 15*(11 - 82*Power(E0,2) + 68*Power(E0,4))*Power(r0,2)))*
      EllipticE(Power(L0,2)/(Power(L0,2) + Power(r0,2))))/(10.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*Pi*Power(r0,10)*Power(Power(L0,2) + Power(r0,2),3.5))
     - (3*(1536*Power(L0,12)*Power(M,2) + 32*Power(L0,10)*M*(593*M - 230*r0)*Power(r0,2) + 
        8*Power(L0,8)*M*Power(r0,4)*(8688*M + (-4255 + 688*Power(E0,2))*r0) + 
        4*Power(L0,6)*Power(r0,6)*(27659*Power(M,2) + 10*(-1426 + 455*Power(E0,2))*M*r0 - 60*(3 - 7*Power(E0,2) + 4*Power(E0,4))*Power(r0,2)) + 
        Power(r0,12)*(680*Power(M,2) + 4*(215 + 227*Power(E0,2))*M*r0 - 15*(47 - 136*Power(E0,2) + 91*Power(E0,4))*Power(r0,2)) + 
        Power(L0,4)*Power(r0,8)*(81664*Power(M,2) + 4*(-10135 + 5139*Power(E0,2))*M*r0 - 15*(143 - 368*Power(E0,2) + 224*Power(E0,4))*Power(r0,2)) + 
        Power(L0,2)*Power(r0,10)*(23772*Power(M,2) + 8*(-1165 + 1096*Power(E0,2))*M*r0 - 15*(142 - 392*Power(E0,2) + 259*Power(E0,4))*Power(r0,2)))*
      EllipticK(Power(L0,2)/(Power(L0,2) + Power(r0,2))))/(20.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*Pi*Power(r0,8)*Power(Power(L0,2) + Power(r0,2),3.5));
}

double Anna_and_Barrys_D_r_circular_orbits(double l, double r1)
{
	return 4.0*(-((M*Power((-2*M + r1)/(-3*M + r1),1.5)*(-34*M + 15*r1)*EllipticE(M/(-2*M + r1)))/
      ((-3 + 4*l + 4*Power(l,2))*Pi*Power(r1,4))) - 
   ((2*M - r1)*Sqrt((-3*M + r1)/(-2*M + r1))*(-110*Power(M,2) + 41*M*r1 + 5*Power(r1,2))*
      EllipticE(M/(-2*M + r1)))/(8.*(-3 + 4*l + 4*Power(l,2))*Pi*(3*M - r1)*Power(r1,4)) + 
   (Sqrt((-2*M + r1)/(-3*M + r1))*(-188*Power(M,2) + 75*M*r1 + Power(r1,2))*EllipticK(M/(-2*M + r1)))/
    (4.*(-3 + 4*l + 4*Power(l,2))*Pi*Power(r1,4)) + 
   (3*Power(-3*M + r1,1.5)*(M + r1)*Hypergeometric2F1(0.5,2.5,1,M/(-2*M + r1)))/
    (16.*(-3 + 4*l + 4*Power(l,2))*Power(r1,4)*Sqrt(-2*M + r1)));
}

double Anna_and_Barrys_F_r(double l, double p, double e, double r1, double ur1, double L){

return (3*(1536*Power(L,18)*M*(8*Power(M,2) + 5*Power(r1,2) + 2*M*r1*(-7 + 60*Power(ur1,2))) - 
        256*Power(L,16)*M*Power(r1,2)*(434*Power(M,2) + Power(r1,2)*(-73 + 210*Power(ur1,2)) - 
           M*r1*(71 + 6660*Power(ur1,2))) + 32*Power(L,14)*Power(r1,4)*
         (-28190*Power(M,3) + 30*Power(r1,3)*(-3 + 2*Power(ur1,2)) + Power(M,2)*r1*(16723 + 188484*Power(ur1,2)) + 
           M*Power(r1,2)*(-1134 - 3563*Power(ur1,2) + 5760*Power(ur1,4))) + 
        4*Power(L,12)*Power(r1,6)*(-618432*Power(M,3) + 2*Power(M,2)*r1*(217735 + 1355376*Power(ur1,2)) + 
           120*Power(r1,3)*(-31 + 12*Power(ur1,2) + 8*Power(ur1,4)) + 
           M*Power(r1,2)*(-55687 + 96826*Power(ur1,2) + 266664*Power(ur1,4))) + 
        Power(r1,18)*(-56616*Power(M,3) + 4*Power(M,2)*r1*(20375 + 10452*Power(ur1,2)) + 
           2*M*Power(r1,2)*(-13328 - 12148*Power(ur1,2) + 3425*Power(ur1,4)) + 
           15*Power(r1,3)*(2 - 28*Power(ur1,2) + 29*Power(ur1,4) + 24*Power(ur1,6))) + 
        Power(L,10)*Power(r1,8)*(-3651616*Power(M,3) + 56*Power(M,2)*r1*(52115 + 189288*Power(ur1,2)) + 
           2*M*Power(r1,2)*(-243083 + 851652*Power(ur1,2) + 1280052*Power(ur1,4)) + 
           15*Power(r1,3)*(-2005 + 482*Power(ur1,2) + 1344*Power(ur1,4) + 128*Power(ur1,6))) + 
        2*Power(L,8)*Power(r1,10)*(-1694836*Power(M,3) + 2*Power(M,2)*r1*(799339 + 1332540*Power(ur1,2)) + 
           4*M*Power(r1,2)*(-86610 + 309665*Power(ur1,2) + 405544*Power(ur1,4)) + 
           15*Power(r1,3)*(-973 + 642*Power(ur1,2) + 1726*Power(ur1,4) + 384*Power(ur1,6))) + 
        Power(L,2)*Power(r1,16)*(-350848*Power(M,3) + 8*Power(M,2)*r1*(63835 + 11676*Power(ur1,2)) + 
           2*M*Power(r1,2)*(-84459 - 19628*Power(ur1,2) + 62863*Power(ur1,4)) + 
           15*Power(r1,3)*(43 + 634*Power(ur1,2) + 1563*Power(ur1,4) + 874*Power(ur1,6))) + 
        Power(L,6)*Power(r1,12)*(-2199264*Power(M,3) + 32*Power(M,2)*r1*(79328 + 31995*Power(ur1,2)) + 
           2*M*Power(r1,2)*(-347026 + 825400*Power(ur1,2) + 1126961*Power(ur1,4)) + 
           15*Power(r1,3)*(-846 + 2708*Power(ur1,2) + 5785*Power(ur1,4) + 1978*Power(ur1,6))) + 
        Power(L,4)*Power(r1,14)*(-1060400*Power(M,3) + 32*Power(M,2)*r1*(45073 + 414*Power(ur1,2)) + 
           2*M*Power(r1,2)*(-227074 + 207132*Power(ur1,2) + 403107*Power(ur1,4)) + 
           15*Power(r1,3)*(-64 + 2312*Power(ur1,2) + 4955*Power(ur1,4) + 2444*Power(ur1,6))))*
      EllipticE(Power(L,2)/(Power(L,2) + Power(r1,2))) - 
     3*Power(r1,2)*(768*Power(L,16)*M*(8*Power(M,2) + 5*Power(r1,2) + 2*M*r1*(-7 + 60*Power(ur1,2))) - 
        32*Power(L,14)*Power(r1,2)*(1904*Power(M,3) + 60*Power(r1,3) + M*Power(r1,2)*(-217 + 1560*Power(ur1,2)) - 
           2*Power(M,2)*r1*(379 + 12330*Power(ur1,2))) - 
        6*Power(L,10)*Power(r1,6)*(229416*Power(M,3) + 15*Power(r1,3)*(207 + 116*Power(ur1,2)) - 
           8*Power(M,2)*r1*(26435 + 78666*Power(ur1,2)) + 2*M*Power(r1,2)*(21088 + 9793*Power(ur1,2) - 42764*Power(ur1,4)))
          - 8*Power(L,12)*Power(r1,4)*(60240*Power(M,3) + 120*Power(r1,3)*(10 + 3*Power(ur1,2)) - 
           2*Power(M,2)*r1*(23795 + 153654*Power(ur1,2)) + M*Power(r1,2)*(6335 + 20386*Power(ur1,2) - 12960*Power(ur1,4)))\
         + Power(L,2)*Power(r1,14)*(-546640*Power(M,3) + 4*Power(M,2)*r1*(149680 + 71721*Power(ur1,2)) + 
           8*M*Power(r1,2)*(-20495 - 12366*Power(ur1,2) + 12042*Power(ur1,4)) + 
           15*Power(r1,3)*(42 + 211*Power(ur1,2) + 280*Power(ur1,4) + 22*Power(ur1,6))) + 
        5*Power(L,8)*Power(r1,8)*(-442288*Power(M,3) + 16*Power(M,2)*r1*(27850 + 39573*Power(ur1,2)) + 
           4*M*Power(r1,2)*(-26362 + 5658*Power(ur1,2) + 50379*Power(ur1,4)) + 
           3*Power(r1,3)*(-1130 - 569*Power(ur1,2) + 272*Power(ur1,4) + 64*Power(ur1,6))) - 
        Power(r1,16)*(91536*Power(M,3) - 4*Power(M,2)*r1*(24740 + 17607*Power(ur1,2)) + 
           4*M*Power(r1,2)*(6664 + 7469*Power(ur1,2) - 1105*Power(ur1,4)) + 
           15*Power(r1,3)*(-2 + 28*Power(ur1,2) + 76*Power(ur1,4) + 81*Power(ur1,6))) + 
        Power(L,6)*Power(r1,10)*(-2243136*Power(M,3) + 4*Power(M,2)*r1*(596812 + 397755*Power(ur1,2)) + 
           4*M*Power(r1,2)*(-155090 + 34235*Power(ur1,2) + 243028*Power(ur1,4)) + 
           15*Power(r1,3)*(-416 + 449*Power(ur1,2) + 1168*Power(ur1,4) + 328*Power(ur1,6))) + 
        3*Power(L,4)*Power(r1,12)*(-482816*Power(M,3) + 4*Power(M,2)*r1*(131848 + 54675*Power(ur1,2)) + 
           4*M*Power(r1,2)*(-35808 - 4135*Power(ur1,2) + 38828*Power(ur1,4)) + 
           5*Power(r1,3)*(24 + 753*Power(ur1,2) + 1252*Power(ur1,4) + 495*Power(ur1,6))))*
      EllipticK(Power(L,2)/(Power(L,2) + Power(r1,2))))/
   (40.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*Pi*Power(r1,12)*Power(Power(L,2) + Power(r1,2),5.5));

}

double Anna_and_Barrys_F_t(double l, double p, double e, double r1, double ur1, double L){

	double e1 = Sqrt(((-2 - 2*e + p)*(-2 + 2*e + p))/(p*(-3 - Power(e,2) + p)));

	return (e1*ur1*(6*(92160*Power(L,18)*Power(M,2) + 3840*Power(L,16)*M*(210*M - r1)*Power(r1,2) + 
          16*Power(L,14)*Power(r1,4)*(182132*Power(M,2) + 60*Power(r1,2) + 9*M*r1*(-43 + 640*Power(ur1,2))) + 
          4*Power(L,12)*Power(r1,6)*(1441156*Power(M,2) + 480*Power(r1,2)*(3 + Power(ur1,2)) + 3*M*r1*(1361 + 44444*Power(ur1,2))) + 
          Power(r1,18)*(27424*Power(M,2) + 2*M*r1*(-7929 + 2155*Power(ur1,2)) + 15*Power(r1,2)*(1 - 15*Power(ur1,2) + 12*Power(ur1,4))) + 
          3*Power(L,10)*Power(r1,8)*(2305324*Power(M,2) + 6*M*r1*(1201 + 70954*Power(ur1,2)) + 5*Power(r1,2)*(977 + 768*Power(ur1,2) + 64*Power(ur1,4))) + 
          Power(L,8)*Power(r1,10)*(5260192*Power(M,2) + 2*M*r1*(-39033 + 807488*Power(ur1,2)) + 15*Power(r1,2)*(1281 + 1966*Power(ur1,2) + 384*Power(ur1,4))) + 
          Power(L,2)*Power(r1,16)*(200372*Power(M,2) + 2*M*r1*(-43371 + 40139*Power(ur1,2)) + 15*Power(r1,2)*(-7 + 201*Power(ur1,2) + 437*Power(ur1,4))) + 
          Power(L,6)*Power(r1,12)*(2599208*Power(M,2) + 2*M*r1*(-101382 + 569383*Power(ur1,2)) + 15*Power(r1,2)*(838 + 2499*Power(ur1,2) + 989*Power(ur1,4))) + 
          3*Power(L,4)*Power(r1,14)*(287728*Power(M,2) + 2*M*r1*(-32108 + 72627*Power(ur1,2)) + 5*Power(r1,2)*(206 + 1389*Power(ur1,2) + 1222*Power(ur1,4))))*
        EllipticE(Power(L,2)/(Power(L,2) + Power(r1,2))) - 3*Power(r1,2)*
        (92160*Power(L,16)*Power(M,2) + 1920*Power(L,14)*M*(387*M - 14*r1)*Power(r1,2) + 
          16*Power(L,12)*Power(r1,4)*(150542*Power(M,2) + 60*Power(r1,2) + 9*M*r1*(-1013 + 720*Power(ur1,2))) + 
          12*Power(L,10)*Power(r1,6)*(349948*Power(M,2) + 10*Power(r1,2)*(41 + 16*Power(ur1,2)) + M*r1*(-29995 + 42444*Power(ur1,2))) + 
          3*Power(L,2)*Power(r1,14)*(145140*Power(M,2) + 4*M*r1*(-13868 + 8773*Power(ur1,2)) + 5*Power(r1,2)*(-15 - 18*Power(ur1,2) + 22*Power(ur1,4))) + 
          3*Power(L,8)*Power(r1,8)*(1471344*Power(M,2) + 4*M*r1*(-45508 + 83005*Power(ur1,2)) + 5*Power(r1,2)*(603 + 656*Power(ur1,2) + 64*Power(ur1,4))) + 
          Power(r1,16)*(72308*Power(M,2) + 4*M*r1*(-7929 + 760*Power(ur1,2)) - 15*Power(r1,2)*(-2 + 30*Power(ur1,2) + 81*Power(ur1,4))) + 
          Power(L,6)*Power(r1,10)*(2993188*Power(M,2) + 4*M*r1*(-141111 + 243553*Power(ur1,2)) + 15*Power(r1,2)*(459 + 1098*Power(ur1,2) + 328*Power(ur1,4))) + 
          3*Power(L,4)*Power(r1,12)*(461284*Power(M,2) + 4*M*r1*(-32833 + 40503*Power(ur1,2)) + 5*Power(r1,2)*(103 + 582*Power(ur1,2) + 495*Power(ur1,4))))*
        EllipticK(Power(L,2)/(Power(L,2) + Power(r1,2)))))/
   (40.*(-3 + 2*l)*(-1 + 2*l)*(3 + 2*l)*(5 + 2*l)*Pi*Power(r1,10)*(-2*M + r1)*Power(Power(L,2) + Power(r1,2),5.5));


}
