// ***************************************************************
// This file was created using the bat-project script
// for project NonPert.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "NonPertFit.h"

int main()
{
    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new NonPertFit object
    NonPertFit m("SIYY-g");

    // set precision
    m.SetPrecision(BCEngineMCMC::kQuick);

    BCLog::OutSummary("Test model created");

    // Add Data
    BCDataSet data;
    data.ReadDataFromFileTxt("ptData.txt",13);
    data.Fix(0); // Fix data type
    data.Fix(1); // Fix QMin
    data.Fix(2); // Fix QMax
    data.Fix(3); // Fix Pt
    data.Fix(6); // Fix dumSys
    data.Fix(7); // Fix dumUnCorr
    data.Fix(9); // Fix <y>
    data.Fix(10); // Fix <Q>
    data.Fix(11); // Fix yk
    data.Fix(12); // Fix efficiency

    // Assign data to the model
    m.SetDataSet(&data);

    //////////////////////////////
    // perform your analysis here

    // Normalize the posterior by integrating it over the full parameter space
    // m.Normalize();

    // Write Markov Chain to a ROOT file as a TTree
    // m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // run mode finding; by default using Minuit
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print summary plots
    // m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
    // m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
    // m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
    // m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");

    // print results of the analysis into a text file
    m.PrintSummary();

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
