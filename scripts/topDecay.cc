#include <vector>
#include <exception>
#include <complex>
#include <cmath>
#include <iostream>
#include <map>

#include "TLorentzVector.h"

#include "../src/utility.hh"
#include "../src/vegas.hh"

#define SMALL 1E-8

std::complex<double> I(0,1);

std::vector<std::complex<double>> yKet(const TLorentzVector &p, const int helicity){
    std::vector<std::complex<double>> results(2);
    double pp = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    if(helicity == 2) {
        results[0] = 1;
        results[1] = 0;
    } else if(helicity == -2) {
        results[0] = 0;
        results[1] = 1.0;
    } else if(helicity == 1) {
        if(pp*(pp+p[2]) < SMALL) {
            results[0] = 0;
            results[1] = -1.0;
        } else {
            double den = sqrt(2.0*pp*(pp+p[2]));
            results[0] = (pp+p[2])/den;
            results[1] = p[0]/den+I*p[1]/den;
        }
    } else if(helicity == -1) {
        if(pp*(pp+p[2]) < SMALL) {
            results[0] = 1.0;
            results[1] = 0;
        } else {
            double den = sqrt(2.0*pp*(pp+p[2]));
            results[1] = (pp+p[2])/den;
            results[0] = -p[0]/den+I*p[1]/den;
        }
    } else {
        throw std::runtime_error("Invalid helicity state");
    }

    return results;

}

std::vector<std::complex<double>> yBra(const TLorentzVector &p, const int helicity){
    std::vector<std::complex<double>> results(2);
    results = yKet(p,helicity);
    results[0] = conj(results[0]);
    results[1] = conj(results[1]);

    return results;
}

std::vector<std::complex<double>> ySlash(const TLorentzVector &p, const int helicity){
    std::vector<std::complex<double>> results(4);
    if(helicity == 1) {
        results[0] = p[3]+p[2];
        results[1] = p[0] - I*p[1];
        results[3] = p[3]-p[2];
    } else if(helicity == -1) {
        results[0] = p[3]-p[2];
        results[1] = -p[0] + I*p[1];
        results[3] = p[3]+p[2];
    } else {
        throw std::runtime_error("Invalid helicity state");
    }
    results[2] = conj(results[1]);

    return results;
}

std::complex<double> yBraKet(const TLorentzVector &p1, const int helicity1, const TLorentzVector &p2, const int helicity2){

    std::vector<std::complex<double>> bra = yBra(p1,helicity1);
    std::vector<std::complex<double>> ket = yKet(p2,helicity2);

    std::complex<double> result = bra[0]*ket[0] + bra[1]*ket[1];

    return result;
}

std::complex<double> yBraSKet(const TLorentzVector &p1, const int helicity1, const TLorentzVector &p3, const int helicity3,
                              const TLorentzVector &p2, const int helicity2){

    std::vector<std::complex<double>> bra = yBra(p1,helicity1);
    std::vector<std::complex<double>> ket = yKet(p2,helicity2);
    std::vector<std::complex<double>> slash = ySlash(p3,helicity3);

    std::complex<double> result = bra[0]*slash[0]*ket[0] + bra[0]*slash[1]*ket[1] + bra[1]*slash[2]*ket[0] + bra[1]*slash[3]*ket[1];

    return result;
}

std::complex<double> yBraSSKet(const TLorentzVector &p1, const int helicity1, const TLorentzVector &p3, const int helicity3,
                               const TLorentzVector &p4, const int helicity4, const TLorentzVector &p2, const int helicity2){

    TLorentzVector zero(0,0,0,0);
    return yBraSKet(p1,helicity1,p3,helicity3,zero,-2)*yBraSKet(zero,2,p4,helicity4,p2,helicity2)
         + yBraSKet(p1,helicity1,p3,helicity3,zero,2)*yBraSKet(zero,-2,p4,helicity4,p2,helicity2);

}

std::vector<double> Calfm3(double Smin, double alphas) {
    double CF = 4.0/3.0;
    double fac = alphas*CF/4.0/M_PI;
    double mTop = 172, mTop2 = pow(mTop,2), mW = 80.38, mW2 = pow(mW,2);
    double W = mW2/mTop2;

    std::vector<double> result(2);
    result[1] = fac/W*log(1-W);
    result[1] *= 2.0/mTop;

    double z1 = W/(W-1);
    double z2 = (1.0-W)/(2.0-W);
    double yy = 2.0*(xLi(2,z1)-xLi(2,z2));
    double x = Smin/mTop2;

    result[0] = fac*(-2.5-2.0/3.0*M_PI*M_PI+pow(log(2.0),2)-2.0*log(2.0)+yy+(3.0+2.0*log(2)-1.0/W)*log(1-W)
              - 2*pow(log(1-W),2)-3.5*log(x)-pow(log(x),2)+(4.0-2*W)/(1-W)*log(2-W)-pow(log(2-W),2)+4*log(1-W)*log(x)
              - 1.0/(1-W));

    return result;
}

double mt;
double gt = 1.465;

double GetAmpTree2(const double wt, std::map<int,TLorentzVector> &momentum, const int TopCharge, const int TopPol) {
    std::complex<double> b1n2 = yBraKet(momentum[5],-TopCharge,momentum[12],TopCharge);
    std::complex<double> e2t1 = yBraKet(momentum[11],TopCharge,momentum[6],-1);
    std::complex<double> e2t2 = yBraKet(momentum[11],TopCharge,momentum[6],+1);


    double termc = momentum[6].Vect().Mag();
    std::map<int,double> con;
    con[1] = sqrt(momentum[6][3]+termc);
    con[-1]= sqrt(momentum[6][3]-termc);

    std::map<int,std::complex<double>> tben;
    tben[1] = -2.0*con[TopCharge]*e2t1*b1n2;
    tben[-1] = -2.0*con[-TopCharge]*e2t2*b1n2;

    double mTop = 172, gTop = 1.465;
    double mW = 80.42, mW2=pow(mW,2), gW=2.085;
    double StProp = 1.0;//M_PI/(mt*gt);
    double SwProp = pow(momentum[24].M2()-mW2,2)+pow(mW*gW,2);
    SwProp = 1.0/SwProp;
    double GMU = 1.166389E-5;
    double gWeak2 = 8.0*GMU*mW2/sqrt(2.0);
    double coupl = pow(gWeak2,2)/4.0*8.0*momentum[5].E()*momentum[11].E()*momentum[12].E();

    return norm(tben[TopPol])*SwProp*StProp*coupl*wt; 
    
}

double GetAmpVirt2(const double wt, std::map<int,TLorentzVector> &momentum, const int TopCharge, const bool TopPol) {
    std::complex<double> b1n2 = yBraKet(momentum[5],-TopCharge,momentum[12],TopCharge);
    std::complex<double> e2t1 = yBraKet(momentum[11],TopCharge,momentum[6],-1);
    std::complex<double> n1b2e1 = yBraSKet(momentum[12],-TopCharge,momentum[5],TopCharge,momentum[11],-TopCharge);
    std::complex<double> b1t1 = yBraKet(momentum[5],-TopCharge,momentum[6],-1);
    std::complex<double> e2t2 = yBraKet(momentum[11],TopCharge,momentum[6],1);
    std::complex<double> b1t2 = yBraKet(momentum[5],-TopCharge,momentum[6],1);
    std::complex<double> b2e1 = yBraKet(momentum[5],TopCharge,momentum[11],-TopCharge);
    std::complex<double> n1t1 = yBraKet(momentum[12],-TopCharge,momentum[6],-1);
    std::complex<double> b2t1 = yBraKet(momentum[5],TopCharge,momentum[6],-1);
    std::complex<double> n1t2 = yBraKet(momentum[12],-TopCharge,momentum[6],1);
    std::complex<double> b2t2 = yBraKet(momentum[5],TopCharge,momentum[6],1);

    std::map<int,std::complex<double>[2]> tben1;
    std::map<int,std::complex<double>[2]> tben2;

    double Smin=5, alphas=0.118;

    std::vector<double> FMBW = Calfm3(Smin,alphas);

    double d1L[2], d2L[2], d1R[2], d2R[2];
    d1L[0] = 1.0;
    d1R[0] = 0.0;
    d2L[0] = 0.0;
    d2R[0] = 0.0;
    d1L[1] = FMBW[0];
    d1R[1] = 0.0;
    d2L[1] = 0.0;
    d2R[1] = FMBW[1];

    double termc = momentum[6].Vect().Mag();
    std::map<int,double> con;
    con[1] = sqrt(momentum[6][3]+termc);
    con[-1]= sqrt(momentum[6][3]-termc);

    for(int i = 0; i < 2; i++){
        std::complex<double> cdump1 = -2.0*d1L[i]*con[(TopCharge-1)/2]*e2t2*b1n2;
        std::complex<double> cdump2 = d2R[i]*con[(1-TopCharge)/2]*n1b2e1*b1t2;
        tben1[i][1] = cdump1 + cdump2;

        cdump1 = -2.0*d1L[i]*con[(1-TopCharge)/2]*e2t1*b1n2;
        cdump2 = d2R[i]*con[(TopCharge-1)/2]*n1b2e1*b1t1;
        tben1[i][0] = cdump1 + cdump2;

        cdump1 = -2.0*d1R[i]*con[-(TopCharge-1)/2]*n1t2*b2e1;
        cdump2 = d2R[i]*con[-(1-TopCharge)/2]*n1b2e1*b2t2;
        tben2[i][1] = cdump1 + cdump2;

        cdump1 = -2.0*d1R[i]*con[-(1-TopCharge)/2]*n1t1*b2e1;
        cdump2 = d2L[i]*con[-(TopCharge-1)/2]*n1b2e1*b2t1;
        tben2[i][0] = cdump1 + cdump2;
    }

    std::complex<double> termd1 = conj(tben1[0][TopPol])*tben1[1][TopPol];
    std::complex<double> termd2 = conj(tben2[0][TopPol])*tben2[1][TopPol];
    std::complex<double> termd = termd1+termd2;

    double mTop = 172, gTop = 1.465;
    double mW = 80.38, mW2=pow(mW,2), gW=2.085;
    double StProp = M_PI/(mTop*gTop);
    double SwProp = pow(momentum[24].M2()-mW2,2)+pow(mW*gW,2);
    SwProp = 1.0/SwProp;
    double GMU = 1.166389E-5;
    double gWeak2 = 8.0*GMU*mW2/sqrt(2.0);
    double coupl = gWeak2/4.0*8.0*momentum[5][3]*momentum[11][3]*momentum[12][3];

    return 2.0*real(termd)*SwProp*StProp*coupl*wt;

}

double GetAmpReal2(const double wt, std::map<int,TLorentzVector> &momentum, const int TopCharge, const int TopPol){
    double mTop = 172, gTop = 1.465;
    double mW = 81.38, mW2=pow(mW,2), gW=2.085;
    double SwProp = pow(momentum[24].M2()-mW2,2)+pow(mW*gW,2);
    SwProp = 1.0/SwProp;
    double GMU = 1.166389E-5;
    double gWeak2 = 8.0*GMU*mW2/sqrt(2.0);
    double coupl = gWeak2/4.0*8.0*momentum[5][3]*momentum[11][3]*momentum[12][3];
    std::complex<double> tDec[2][2];
    double StProp = M_PI/(mTop*gTop);
    double SbProp = 1.0/momentum[18].M2();
    double termc = momentum[6].Vect().Mag();
    std::map<int,double> con;
    con[1] = sqrt(momentum[6][3]+termc);
    con[-1]= sqrt(momentum[6][3]-termc);
    for(int i = 0; i < 2; i++) {
        TLorentzVector gluon;
        if(i == 0) gluon = momentum[21];
        if(i == 1) gluon = momentum[9];

        std::complex<double> TM1 = yBraKet(momentum[5],-TopCharge,momentum[12],TopCharge);
        std::complex<double> TM2 = yBraKet(momentum[11],TopCharge,momentum[6],-1);
        std::complex<double> TM3 = yBraKet(momentum[11],TopCharge,momentum[6],+1);

        std::complex<double> TM11 = yBraSSKet(momentum[11],TopCharge,momentum[17],-TopCharge,gluon,TopCharge,momentum[6],1);
        std::complex<double> TM12 = yBraSKet (momentum[11],TopCharge,gluon,-TopCharge,momentum[6],1);
        std::complex<double> TM13 = yBraSSKet(momentum[11],TopCharge,momentum[17],-TopCharge,gluon,TopCharge,momentum[6],-1);
        std::complex<double> TM14 = yBraSKet (momentum[11],TopCharge,gluon,-TopCharge,momentum[6],-1);
        std::complex<double> TM21 = yBraSSKet(momentum[5],-TopCharge,gluon,TopCharge,momentum[18],-TopCharge,momentum[12],TopCharge);

        std::complex<double> TBeng1[2], TBeng2[2];
        TBeng1[1] = 2.0*TM1*(con[(TopCharge-1)/2]*TM11+mTop*con[(1-TopCharge)/2]*TM12);
        TBeng2[1] = 2.0*con[(TopCharge-1)/2]*TM2*TM21;
        TBeng1[0] = 2.0*TM1*(con[(1-TopCharge)/2]*TM13+mTop*con[(TopCharge-1)/2]*TM14);
        TBeng2[0] = 2.0*con[(1-TopCharge)/2]*TM3*TM21;

        tDec[i][0] = TBeng1[0]*StProp+TBeng2[0]*SbProp;
        tDec[i][1] = TBeng1[1]*StProp+TBeng2[1]*SbProp;
    }

    double amp2 = 0;
    for(int i = 0; i < 2; i++) {
        amp2 += norm(tDec[i][TopPol]);
    }

    return amp2*SwProp*StProp*coupl*wt;
    
}

double x2Phase(double *x,double m12,double m1,double m2,TLorentzVector &p1,TLorentzVector &p2) {
    double costh=2.0*x[0]-1.0;
    double sinth=sqrt(1-costh*costh);
    double phiSta=2.0*M_PI*x[1];
    double rsh=m12;
    double sh=rsh*rsh;
    double temp1=(sh+m1*m1-m2*m2)/2.0/rsh;
    double q1cm=sqrt(temp1*temp1-m1*m1);
    double sphi1=sin(phiSta);
    double cphi1=cos(phiSta);
    p1.SetPxPyPzE(q1cm*cphi1*sinth,q1cm*sphi1*sinth,q1cm*costh,temp1);
    p2=-p1;
    p2[3]=m12-p1[3];

    double dPhase2=q1cm/(4.0*rsh*pow(2.0*M_PI,6));
    double yJaco1=2.0;
    double yJaco2=2.0*M_PI;
    return yJaco1*yJaco2*pow(2.0*M_PI,4)*dPhase2;
}

double calcDecay(double *x, double wgt) {
    std::map<int,TLorentzVector> momentum;
    momentum[6].SetPxPyPzE(0,0,0,mt);
    double m1 = 0;

    double mW = 80.42, mW2 = pow(mW,2), gW = 2.085;
    double yMwMin = (-mW2)/(mW*gW);
    double yMwMax = (pow(mt,2)-mW2)/(mW*gW);
    double yMwVar = atan(yMwMax)*x[4]+(1-x[4])*atan(yMwMin);
    double yMwNorm = atan(yMwMax)-atan(yMwMin);
    double m22 = mW2 + mW*gW*atan(yMwVar);
    double m2 = sqrt(m22);
    double wt0 = mW*gW*yMwNorm/pow(cos(yMwVar),2);

    double x1[2] = {x[0],x[1]};
    double wt1 = x2Phase(x1,momentum[6].M(),m1,m2,momentum[5],momentum[24]);
    momentum[5].Boost(momentum[6].BoostVector());
    momentum[24].Boost(momentum[6].BoostVector());
    double x2[2] = {x[2],x[3]};
    double wt2 = x2Phase(x2,momentum[24].M(),0,0,momentum[11],momentum[12]);
    momentum[11].Boost(momentum[24].BoostVector());
    momentum[12].Boost(momentum[24].BoostVector());

    double spinAvg = 0.5;
    double flux = 2.0*mt;
    double GMU = 1.166389E-5;
    double gWeak2 = 8.0*GMU*mW2/sqrt(2.0);
    double WDecay = gWeak2*mW/(48*M_PI);

    return wt0*wt1*wt2*(GetAmpTree2(1,momentum,1,1)+GetAmpTree2(1,momentum,1,-1))*spinAvg/flux/(2*M_PI)/WDecay;
}

int main(int argc, char** argv) {

    mt = atof(argv[1]);

    std::map<int,TLorentzVector> momentum;
    momentum[6].SetPxPyPzE(0,0,0,mt);
    momentum[5].SetPxPyPzE(0.2914595E2,0.6530043E2,-0.3328243E2,0.7887555E2);
    momentum[11].SetPxPyPzE(-0.3152738E1,-0.5980514E2,0.1829213E2,0.6261946E2);
    momentum[12].SetPxPyPzE(-0.2599322E2,-0.5495290E1,0.1499029E2,0.3050499E2);
    momentum[24] = momentum[11]+momentum[12];
    momentum[0] = momentum[24]+momentum[5];
    double mW = 80.42, mW2=pow(mW,2), gW=2.085;
    double GMU = 1.166389E-5;
    double gWeak2 = 8.0*GMU*mW2/sqrt(2.0);
    double WDecay = gWeak2*mW/(48*M_PI);
    std::cout << WDecay << std::endl;
    std::cout << momentum[5].Px() << "\t" << momentum[5].Py() << "\t" << momentum[5].Pz() << "\t" << momentum[5].E() << std::endl;
    std::cout << momentum[11].Px() << "\t" << momentum[11].Py() << "\t" << momentum[11].Pz() << "\t" << momentum[11].E() << std::endl;
    std::cout << momentum[12].Px() << "\t" << momentum[12].Py() << "\t" << momentum[12].Pz() << "\t" << momentum[12].E() << std::endl;
    std::cout << momentum[24].Px() << "\t" << momentum[24].Py() << "\t" << momentum[24].Pz() << "\t" << momentum[24].E() << "\t" << momentum[24].M() << std::endl;

    std::cout << GetAmpTree2(1,momentum,1,1)  << std::endl;
    std::cout << GetAmpTree2(1,momentum,1,-1)  << std::endl;
    std::cout << GetAmpTree2(1,momentum,-1,-1) << std::endl;
    std::cout << GetAmpTree2(1,momentum,-1,1) << std::endl;

    double a[5] = {0,0,0,0,0}, b[5] = {1,1,1,1,1};
    Vegas v(5,46648643,20,2000000);
    v.nprn=0;
    v.SetRange(a,b);
    double intGT = v.VEGAS(calcDecay);
    v.init=1;
    v.itmx=40;
    v.SetNCalls(4000000);
    intGT = v.VEGAS(calcDecay);
    double mTop = 172, gTop = 1.465;
    double beta = sqrt(1-mW2/pow(mt,2));
    double GamT = GMU*pow(mt,3)/(8*sqrt(2)*M_PI)*pow(beta,4)*(3-2*pow(beta,2));
    std::cout << intGT << "\t" << GamT << "\t" << intGT/GamT << std::endl;

    return 0;
}
