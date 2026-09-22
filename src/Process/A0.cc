#include <cmath>

#include "ResBos/Beam.hh"
#include "ResBos/Electroweak.hh"
#include "ResBos/PhaseSpace.hh"
#include "ResBos/Process/DrellYan.hh"
#include "ResBos/QCDConst.hh"
#include "ResBos/ResBos.hh"
#include "ResBos/Settings.hh"

#include "ResBos/loguru.hpp"

namespace ResBos {

REGISTER(A0)

A0::A0(std::shared_ptr<ResBos> resbos_, const IO::Settings &settings) : Process(resbos_,settings) {
    name = "A0";
}

// Process specific initialization
void A0::Initialize([[maybe_unused]]const IO::Settings &settings) {
#ifndef FITTING
    ps = std::unique_ptr<PhaseSpace>(new PhaseSpace(pow(resbos -> GetECM(),2)));
    std::vector<double> Mass = {0,0,0};
    ps->SetMass(Mass);
    ps->SetQRange(settings.GetSettingDouble("Qmin"), settings.GetSettingDouble("Qmax"));
    ps->SetQtRange(settings.GetSettingDouble("Qtmin"), settings.GetSettingDouble("Qtmax"));
#endif // FITTING

    double FAC_A = sqrt(2.0*abs(GetEW() -> GetACoupl()))/2.0;

    double QL_A1 = (2.0/3.0)*FAC_A;
    double QR_A1 = (2.0/3.0)*FAC_A;
    double QL_A2 = (-1.0/3.0)*FAC_A;
    double QR_A2 = (-1.0/3.0)*FAC_A;
    double EL_A = -1.0*FAC_A;
    double ER_A = -1.0*FAC_A;

    // 0 corresponds to up-type quarks, 1 corresponds to down-type quarks
    C_E[0] = ( QL_A1*QL_A1 + QR_A1*QR_A1 )*( EL_A*EL_A + ER_A*ER_A ); 
    C_O[0] = ( QL_A1*QL_A1 - QR_A1*QR_A1 )*( EL_A*EL_A - ER_A*ER_A );
    C_E[1] = ( QL_A2*QL_A2 + QR_A2*QR_A2 )*( EL_A*EL_A + ER_A*ER_A ); 
    C_O[1] = ( QL_A2*QL_A2 - QR_A2*QR_A2 )*( EL_A*EL_A - ER_A*ER_A );
    DLOG_F(INFO,"Photon Coupling to up quarks and electrons: %f,%f",C_E[0],C_O[0]);
    DLOG_F(INFO,"Photon Coupling to down quarks and electrons: %f,%f",C_E[1],C_O[1]);
}

double A0::Charge(const int &pid, const size_t&) const {
    return charge.at(abs(pid));
}

double A0::ME(const std::vector<double> &psPoint, const std::vector<double> &wgt) const {
    // psPoint[0] = Q
    // psPoint[1] = q_T
    // psPoint[2] = y
    // psPoint[3] = theta_star
    // psPoint[4] = phi_star (Not used in resummed or asymptotic piece)
    const double EW_ALFA = 1.0/137.0359895;

    double AA_PROP = 1.0;

    double ALFA_Q2 = 1.0/132.3572536;
    AA_PROP = AA_PROP*pow(ALFA_Q2/EW_ALFA,2);

    double SigS[2] = {wgt[0],wgt[2]};
    double SigA[2] = {wgt[1],wgt[3]};
    double the_sta = psPoint[3];

    double ANGFUNC_ODD[2], ANGFUNC_EVEN[2];
    ANGFUNC_EVEN[0] = (1.0 + pow(cos(the_sta),2))*(SigS[0]);
    ANGFUNC_ODD[0] = (2.0*cos(the_sta))*(SigA[0]);
    ANGFUNC_EVEN[1] = (1.0 + pow(cos(the_sta),2))*(SigS[1]);
    ANGFUNC_ODD[1] = (2.0*cos(the_sta))*(SigA[1]);

    double COEFF_EVEN[2];
    COEFF_EVEN[0] = AA_PROP*C_E[0];
    COEFF_EVEN[1] = AA_PROP*C_E[1];
    double COEFF_ODD[2];
    COEFF_ODD[0] = AA_PROP*C_O[0];
    COEFF_ODD[1] = AA_PROP*C_O[1];

    double FACTOR=2.0*(M_PI/3.0)*(4.0*2.0);

    return FACTOR*(( COEFF_EVEN[0]*ANGFUNC_EVEN[0] + COEFF_ODD[0]*ANGFUNC_ODD[0])+
            ( COEFF_EVEN[1]*ANGFUNC_EVEN[1] + COEFF_ODD[1]*ANGFUNC_ODD[1]));
}

double A0::GetCoupling(const double &Q) const {
    return M_PI*pow(alpha(Q),2)*4.0/9.0/Q/Q;
}

Partons A0::GetPartons(const size_t&) const {
    Partons partons;
    partons.push_back(std::make_pair(2,-2));
    partons.push_back(std::make_pair(4,-4));
    partons.push_back(std::make_pair(1,-1));
    partons.push_back(std::make_pair(3,-3));
    partons.push_back(std::make_pair(5,-5));
    return partons;
}

double A0::CxFCxF(const double &mu, const double &x1, const double &x2, const size_t &mode) const {
    double total1 = 0, total2 = 0;
    Partons partons = GetPartons(mode/2);
    for(auto parton : partons) {
        total1 += charge.at(parton.first)
               * resbos -> GetBeams().first -> PDF(Beam::BDF::C,parton.first,x1,mu)
               * resbos -> GetBeams().second -> PDF(Beam::BDF::C,parton.second,x2,mu);
        total2 += charge.at(parton.first)
               * resbos -> GetBeams().first -> PDF(Beam::BDF::C,parton.second,x1,mu)
               * resbos -> GetBeams().second -> PDF(Beam::BDF::C,parton.first,x2,mu);
    }

    if(mode % 2 == 0) return total1+total2;
    return total1 - total2;
}

double A0::GetPhaseSpacePoint(const std::vector<double> &x, std::vector<double>& Point,
                              std::vector<FourVector>& momenta) const {
    return ps -> Wresph(x, 0.0, Point, momenta);
}

double A0::H1() const { 
    return QCD::CF*(QCD::pi2/2.0-4); 
}

double A0::H2(const double &Q) const {
    const size_t nf = resbos -> GetNF(Q);
    return H2(nf);
}

double A0::H2(const size_t &nf) const {
    return QCD::CF*QCD::CA*(59.0*QCD::ZETA3/18.0-1535.0/192.0+215.0*QCD::pi2/216.0
            -pow(QCD::pi2,2)/240.0)
        + 0.25*QCD::CF*QCD::CF*(-15*QCD::ZETA3+511.0/16.0-67*QCD::pi2/12.0
                +17*pow(QCD::pi2,2)/45.0)
        + QCD::CF*static_cast<double>(nf)/864.0*(192*QCD::ZETA3+1143-152*QCD::pi2);
}

double A0::GetPert(const double &x1, const double &x2, const double &sh, const double &th,
                   const double &uh, const double &mu, const double &Q,
                   const double &qt, const size_t &mode) const {
    double qqb = 0, qg = 0, gq = 0, qbg = 0, gqb = 0;
    double gPDF1 = resbos -> GetBeams().first -> PDF(Beam::BDF::PDF,0,x1,mu);
    double gPDF2 = resbos -> GetBeams().second -> PDF(Beam::BDF::PDF,0,x2,mu);
    int signQQB = 1, signG1 = 1, signG2 = 1, signG3 = 1, signG4 = 1;
    if(mode % nAngular == 1) {
        signG2 = -1;
        signG4 = -1;
    } else if(mode % nAngular == 3) {
        signQQB = -1;
        signG1 = -1;
        signG4 = -1;
    } else if(mode % nAngular == 4) {
        signG1 = -1;
        signG2 = -1;
    }
    Partons partons = GetPartons(mode/nAngular);
    for(auto parton : partons) {
        qqb +=  1.0/sh*(
                charge.at(parton.first) *
                resbos -> GetBeams().first -> PDF(Beam::BDF::PDF,parton.first,x1,mu)*
                resbos -> GetBeams().second -> PDF(Beam::BDF::PDF,parton.second,x2,mu)*
                Amp2QQ(sh,th,uh,Q,qt,mode) + signQQB*
                charge.at(parton.first) *
                resbos -> GetBeams().first -> PDF(Beam::BDF::PDF,parton.second,x1,mu)*
                resbos -> GetBeams().second -> PDF(Beam::BDF::PDF,parton.first,x2,mu)*
                Amp2QQ(sh,uh,th,Q,qt,mode));

        gqb += gPDF1*resbos -> GetBeams().second -> PDF(Beam::BDF::PDF,parton.second,x2,mu)
             * charge.at(parton.first);
        gq  += gPDF1*resbos -> GetBeams().second -> PDF(Beam::BDF::PDF,parton.first,x2,mu)
             * charge.at(parton.first);

        qg  += gPDF2*resbos -> GetBeams().first -> PDF(Beam::BDF::PDF,parton.first,x1,mu)
             * charge.at(parton.first);
        qbg += gPDF2*resbos -> GetBeams().first -> PDF(Beam::BDF::PDF,parton.second,x1,mu)
             * charge.at(parton.first);
    }

    double gContrib = 1.0/sh * (
            (signG1*gq + signG3*gqb) * Amp2GQ(sh,th,uh,Q,qt,mode) +
            (signG2*qg + signG4*qbg) * Amp2GQ(sh,uh,th,Q,qt,mode));

    return qqb + gContrib;
}

double A0::Amp2QQ(const double &sh, const double &th, const double &uh,
                  const double &Q, const double &qt, const size_t &mode) const {
    double result = (pow(sh+uh,2) + pow(sh+th,2))/uh/th; // = t/u + u/t + 2sQ^2/(ut)
    double mt = sqrt(Q*Q + qt*qt);

    double factor = 1.0;
    if(mode % nAngular == 1) 
        factor = qt*Q/mt/mt*(pow(Q*Q-uh,2)-pow(Q*Q-th,2))/(pow(Q*Q-uh,2)+pow(Q*Q-th,2));
    else if(mode % nAngular == 2)
        factor = qt*qt/mt/mt;
    else if(mode % nAngular == 3)
        factor = Q/mt;
    else if(mode % nAngular == 4)
        factor = 2*qt/mt*(pow(Q*Q-uh,2)-pow(Q*Q-th,2))/(pow(Q*Q-uh,2)+pow(Q*Q-th,2));

    return factor*result*2.0/3.0/M_PI;
}

double A0::Amp2GQ(const double &sh, const double &th, const double &uh,
                  const double &Q, const double &qt, const size_t &mode) const {
    double result = -(pow(sh+th,2)+pow(uh+th,2))/uh/sh/4.0/M_PI;
    double mt = sqrt(Q*Q + qt*qt);

    double factor = 1.0;
    if(mode % nAngular == 1) 
        factor = qt*Q/mt/mt*(2.0*pow(sh+th,2)-pow(sh+uh,2))/(pow(Q*Q-sh,2)+pow(Q*Q-uh,2));
    else if(mode % nAngular == 2)
        factor = qt*qt/mt/mt*(pow(Q*Q+sh,2)+pow(Q*Q-uh,2))/(pow(Q*Q-sh,2)+pow(Q*Q-uh,2));
    else if(mode % nAngular == 3)
        factor = Q/mt*(1.0 - 2.0*uh*(Q*Q-sh)/(pow(Q*Q-sh,2)+pow(Q*Q-uh,2)));
    else if(mode % nAngular == 4)
        factor = 2*qt/mt*(1.0+(2.0*sh*(Q*Q-sh))/(pow(Q*Q-sh,2)+pow(Q*Q-uh,2)));

    return factor*result;
}

void A0::MakeEvent(Event &event, const std::vector<FourVector>& momenta) {
    Particle ABoson(pids[0],momenta[0],2), 
             lepton(pids[1],momenta[1],1),
             alepton(pids[2],momenta[2],1);

    ABoson.SetDaughters(2,3);
    lepton.SetMothers(1,0);
    alepton.SetMothers(1,0);

    event.AddParticle(ABoson);
    event.AddParticle(lepton);
    event.AddParticle(alepton);
}

double A0::alpha(const double &Q) const {
    return GetEW() -> GetAlphaEM(Q*Q);
}


}
