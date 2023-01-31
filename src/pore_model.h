#ifndef POREMODEL_H_
#define POREMODEL_H_

#include <string>
#include <vector>

namespace sigmap {
struct PoreModelParameters {
  uint16_t kmer;
  float level_mean;
  float level_stdv;
  float sd_mean;
  float sd_stdv;
  // float weight; // No need to store weight
};

class PoreModel {
 public:
  PoreModel() {}
  ~PoreModel() {}
  void Load(const std::string &pore_model_file_path);
  void Print();
  std::vector<float> GetLevelMeansAt(const char *sequence,
                                     uint32_t start_position,
                                     uint32_t end_position) const;
  int GetKmerSize() const { return kmer_size_; }

  std::vector<PoreModelParameters> pore_models_;

 protected:
  int kmer_size_;
};
}  // namespace sigmap

#endif  // POREMODEL_H_
