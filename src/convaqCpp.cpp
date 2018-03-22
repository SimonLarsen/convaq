#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <random>
#include <thread>
#include "defines.h"
#include "Segment.h"
#include "Region.h"
#include "CNVR.h"
#include "df_to_segments.h"
#include "get_regions.h"
#include "statistical_model.h"

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame convaqCpp(
    DataFrame df1,
    DataFrame df2,
    uint model_num,
    double cutoff,
    bool qvalues,
    uint qvalues_rep,
    uint nthreads
) {
  if(nthreads == 0) nthreads = std::thread::hardware_concurrency();

  MODEL model = (MODEL)model_num;

  std::vector<Segment> segments1, segments2;

  // convert data frames to vector of Segment objects
  df_to_segments(df1, segments1);
  df_to_segments(df2, segments2);

  // get number of patients in each group
  IntegerVector patients1 = df1["patient"];
  IntegerVector patients2 = df2["patient"];
  int npatients1 = Rcpp::max(patients1)+1;
  int npatients2 = Rcpp::max(patients2)+1;

  // grab all chromosomes in data set
  std::unordered_set<std::string> chromosomes;
  StringVector chr1 = df1["chr"];
  StringVector chr2 = df2["chr"];
  for(size_t i = 0; i < chr1.length(); ++i) chromosomes.insert(std::string(chr1[i]));
  for(size_t i = 0; i < chr2.length(); ++i) chromosomes.insert(std::string(chr2[i]));

  std::vector<Region> regions;
  get_regions(segments1, segments2, npatients1, npatients2, chromosomes, regions);

  std::vector<CNVR> results;

  if(model == MODEL_STAT) {
    statistical_model(regions, npatients1, npatients2, cutoff, results);
  }

  // sort by p-value
  std::sort(results.begin(), results.end(), [](CNVR &a, CNVR &b) { return a.pvalue < b.pvalue; });

  if(qvalues) {
    std::vector<std::vector<int>> best(4);

    for(size_t i = 0; i < 4; ++i) best[i].resize(qvalues_rep, 0);

    std::vector<std::thread> threads;
    for(size_t tid = 0; tid < nthreads; ++tid) {
      threads.push_back(std::thread([&best,&segments1,&segments2,&chromosomes,model,cutoff,qvalues_rep,npatients1,npatients2,nthreads](size_t offset) {
        // collect all group-patient pairs
        std::vector<std::pair<int,int>> all_patients;
        for(int i = 0; i < npatients1; ++i) all_patients.emplace_back(0, i);
        for(int i = 0; i < npatients2; ++i) all_patients.emplace_back(1, i);

        std::random_device rd;
        std::minstd_rand rand(rd());

        for(size_t rep = offset; rep < qvalues_rep; rep += nthreads) {
          std::shuffle(all_patients.begin(), all_patients.end(), rand);

          std::vector<std::set<int>> selected(2);
          for(size_t i = 0; i < npatients1; ++i) {
            selected[all_patients[i].first].insert(all_patients[i].second);
          }

          std::vector<Segment> q_segments1, q_segments2;
          for(const Segment &s : segments1) {
            if(selected[0].find(s.patient) != selected[0].end()) {
              q_segments1.push_back(s);
            } else {
              q_segments2.push_back(s);
            }
          }
          for(const Segment &s : segments2) {
            if(selected[1].find(s.patient) != selected[1].end()) {
              q_segments1.push_back(s);
            } else {
              q_segments2.push_back(s);
            }
          }

          std::vector<Region> q_regions;
          get_regions(q_segments1, q_segments2, npatients1, npatients2, chromosomes, q_regions);

          std::vector<CNVR> q_results;
          if(model == MODEL_STAT) {
            statistical_model(q_regions, npatients1, npatients2, cutoff, q_results);
          } else if(model == MODEL_QUERY) {

          }

          for(const CNVR &r : q_results) {
            best[r.type][rep] = std::max(best[r.type][rep], r.length);
          }
        }
      }, tid));
    }

    for(std::thread &th : threads) th.join();

    for(CNVR &r : results) {
      int better = 0;
      for(int l : best[r.type]) {
        if(l > r.length) ++better;
      }
      r.qvalue = (double)better / qvalues_rep;
    }
  }

  // convert CNVR list to data frame
  size_t n = results.size();

  std::vector<std::string> df_chr;
  std::vector<int> df_start, df_end, df_length, df_type;
  std::vector<double> df_pvalue, df_qvalue;

  std::transform(results.begin(), results.end(), std::back_inserter(df_chr),    [](const CNVR &r){ return r.chr; });
  std::transform(results.begin(), results.end(), std::back_inserter(df_start),  [](const CNVR &r){ return r.start; });
  std::transform(results.begin(), results.end(), std::back_inserter(df_end),    [](const CNVR &r){ return r.end; });
  std::transform(results.begin(), results.end(), std::back_inserter(df_length), [](const CNVR &r){ return r.length; });
  std::transform(results.begin(), results.end(), std::back_inserter(df_type),   [](const CNVR &r){ return r.type; });
  std::transform(results.begin(), results.end(), std::back_inserter(df_pvalue), [](const CNVR &r){ return r.pvalue; });
  std::transform(results.begin(), results.end(), std::back_inserter(df_qvalue), [](const CNVR &r){ return r.qvalue; });

  return DataFrame::create(
    Named("chr")    = df_chr,
    Named("start")  = df_start,
    Named("end")    = df_end,
    Named("length") = df_length,
    Named("type")   = df_type,
    Named("pvalue") = df_pvalue,
    Named("qvalue") = df_qvalue
  );
}
