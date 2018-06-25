#include <stdlib.h>
#include <TStopwatch.h>
#include "UnigenQA.h"

using namespace std;

void RunModelQA (TString filePath = "/home/ogolosov/Desktop/analysis/mc/root/dcmqgsm_12.root", 
								 double plab = 12., 
								 TString qaPath = "/home/ogolosov/Desktop/analysis/mc/root/qa_botvina_T0_12agev.root") 
{
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	
	qa::UnigenQA qa;
	qa.Init (filePath, "events");
	qa.SetPlab (plab);
	qa.Init_Histograms();
	qa.Run();
	qa.Write_Histograms(qaPath);    
	
	timer.Stop();
	printf("Real time: %f\n",timer.RealTime());
	printf("CPU time: %f\n",timer.CpuTime());
}

# ifndef __CINT__
int main (int argc, char **argv) {
		TString in = argv [1];
		double pbeam = atof (argv [2]);
		TString out = argv [3];
		
		in = "/home/ogolosov/Desktop/analysis/mc/root/dcmqgsm_12.root";
		pbeam = 13.; 
		out = "/home/ogolosov/Desktop/analysis/mc/qa/qa_botvina_T0_12agev.root";
		
    RunModelQA (in, pbeam, out);
    return 0;
}
# endif