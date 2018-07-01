#include <stdlib.h>
#include <string>
#include <TStopwatch.h>
#include "UnigenQA.h"

using namespace std;

void RunModelQA(TString filePath = "/home/ogolosov/Desktop/analysis/mc/root/dcmqgsm_12.root",
                TString qaPath = "/home/ogolosov/Desktop/analysis/mc/root/qa_botvina_T0_12agev.root") {
  TStopwatch timer;
  timer.Reset();
  timer.Start();

  qa::UnigenQA qa;
  qa.Init(filePath, "events");
  qa.Init_Histograms();
  qa.Run();
  qa.Write_Histograms(qaPath);

  timer.Stop();
  printf("Real time: %f\n", timer.RealTime());
  printf("CPU time: %f\n", timer.CpuTime());
}

# ifndef __CINT__

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv) {

  po::options_description ops_desc("Allowed options");
  ops_desc.add_options()
      ("input,i", po::value<std::string>()->required(), "Input file")
      ("output,o", po::value<std::string>()->required(), "Output file");

  po::positional_options_description pos_desc;
  pos_desc
      .add("input", 1)
      .add("output", 1);

  po::variables_map var_map;
  po::store(po::command_line_parser(argc, argv).options(ops_desc).positional(pos_desc).run(), var_map);
  po::notify(var_map);

  TString in;
  if (var_map.count("input")) {
    string input = var_map["input"].as<string>();
    in = input.c_str();
  }

  TString out;
  if (var_map.count("output")) {
    string output = var_map["output"].as<string>();
    out = output.c_str();
  }

  RunModelQA(in, out);
  return 0;
}
# endif
