#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <random>
#include <thread>
#include <boost/format.hpp>
#include "defines.h"
#include "Segment.h"
#include "Region.h"
#include "CNVR.h"
#include "df_to_segments.h"
#include "get_regions.h"
#include "statistical_model.h"
#include "query_model.h"

using namespace Rcpp;

// [[Rcpp::export]]
List convaqCpp(
    DataFrame df1,
    DataFrame df2,
    uint model_num,
    bool qvalues,
    uint qvalues_rep,
    uint nthreads,
    double cutoff,
    uint comp1, double value1, uint eq1, uint type1,
    uint comp2, double value2, uint eq2, uint type2
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
  } else {
    query_model(
      regions, npatients1, npatients2,
      (COMPARISON)comp1, value1, (EQUALITY)eq1, (VARIATION_TYPE)type1,
      (COMPARISON)comp2, value2, (EQUALITY)eq2, (VARIATION_TYPE)type2,
      results
    );
  }
  
  if(results.size() == 0) {
    return NULL;
  }

  // sort by p-value
  std::sort(results.begin(), results.end(), [](const CNVR &a, const CNVR &b) { return a.pvalue < b.pvalue; });

  if(qvalues) {
    std::vector<std::vector<int>> best(4);

    for(size_t i = 0; i < 4; ++i) best[i].resize(qvalues_rep, 0);

    std::vector<std::thread> threads;
    for(size_t tid = 0; tid < nthreads; ++tid) {
      threads.push_back(std::thread([&](size_t offset) {
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
            query_model(
              q_regions, npatients1, npatients2,
              (COMPARISON)comp1, value1, (EQUALITY)eq1, (VARIATION_TYPE)type1,
              (COMPARISON)comp2, value2, (EQUALITY)eq2, (VARIATION_TYPE)type2,
              q_results
            );
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
    
    std::sort(results.begin(), results.end(), [](const CNVR &a, const CNVR &b) { return a.qvalue < b.qvalue; });
  }

  // prepare output
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
  
  List out = List::create(
    Named("chr") = df_chr,
    Named("start") = df_start,
    Named("end") = df_end,
    Named("length") = df_length,
    Named("type") = df_type,
    Named("pvalue") = df_pvalue,
    Named("qvalue") = df_qvalue
  );
  
  // add within group type frequencies
  std::vector<double> freq(results.size());
  for(size_t group = 0; group < 2; ++group) {
    for(size_t type = 0; type < 3; ++type) {
      for(size_t i = 0; i < results.size(); ++i) {
        freq[i] = regions[results[i].region].get_freq(group, type);
      }
      out.push_back(freq, (boost::format("freq_%d_%d") % group % type).str());
    }
  }

  return out;
}
