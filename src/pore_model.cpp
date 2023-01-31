#include "pore_model.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>

// For sequence manipulation
static constexpr uint8_t char_to_uint8_table_[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
static constexpr char uint8_to_char_table_[8] = {'A', 'C', 'G', 'T', 'N', 'N', 'N', 'N'};

uint8_t CharToUint8(const char c) {
  return char_to_uint8_table_[(uint8_t)c];
}

char Uint8ToChar(const uint8_t i) {
  return uint8_to_char_table_[i];
}

uint64_t GenerateSeedFromSequence(const char *sequence, size_t sequence_length, uint32_t start_position, uint32_t seed_length) {
  uint64_t mask = (((uint64_t)1) << (2 * seed_length)) - 1;
  uint64_t seed = 0;
  for (uint32_t i = 0; i < seed_length; ++i) {
    if (start_position + i < sequence_length) {
      uint8_t current_base = CharToUint8(sequence[i + start_position]);
      if (current_base < 4) {                        // not an ambiguous base
        seed = ((seed << 2) | current_base) & mask;  // forward k-mer
      } else {
        seed = (seed << 2) & mask;  // N->A
      }
    } else {
      seed = (seed << 2) & mask;  // Pad A
    }
  }
  return seed;
}

namespace sigmap {
void PoreModel::Load(const std::string &pore_model_file_path) {
  
//   double real_start_time = ri_realtime();
  int num_kmers = 0;
  std::ifstream pore_model_file_stream(pore_model_file_path);
  std::string pore_model_file_line;
  bool first_line = true;
  while (getline(pore_model_file_stream, pore_model_file_line)) {
    std::stringstream pore_model_file_line_string_stream(pore_model_file_line);
    // skip the header
    if (pore_model_file_line[0] == '#' ||
        pore_model_file_line.find("kmer") == 0) {
      continue;
    }
    std::string kmer;
    pore_model_file_line_string_stream >> kmer;
    if (first_line) {
      kmer_size_ = kmer.length();
      // Allocate memory to save pore model parameters
      size_t num_pore_models = 1 << (kmer_size_ * 2);
      pore_models_.assign(num_pore_models, PoreModelParameters());
      first_line = false;
    }
    assert(kmer.length() == (size_t)kmer_size_);
	//canf: These are 2-bit encoded values of k-mers. K-mers are 5/6-mer long so 2^(2*k) many possible values.
	//canf: pore_models_ hold the pre-determeined kmer models.
    uint64_t kmer_hash_value = GenerateSeedFromSequence(kmer.data(), kmer_size_, 0, kmer_size_);
    PoreModelParameters &pore_model_parameters = pore_models_[kmer_hash_value];
    pore_model_file_line_string_stream >> pore_model_parameters.level_mean >>
        pore_model_parameters.level_stdv >> pore_model_parameters.sd_mean >>
        pore_model_parameters.sd_stdv;
    ++num_kmers;
  }
  std::cerr << "Loaded " << num_kmers << " kmers from the pore model file.\n";
}

void PoreModel::Print() {
  size_t num_pore_models = 1 << (kmer_size_ * 2);
  for (size_t i = 0; i < num_pore_models; ++i) {
    std::cerr << i << " " << pore_models_[i].level_mean << " "
              << pore_models_[i].level_stdv << " " << pore_models_[i].sd_mean
              << " " << pore_models_[i].sd_stdv << "\n";
  }
}

// This funtion return a array of level means given the [start, end) positions
// (0-based).
std::vector<float> PoreModel::GetLevelMeansAt(const char *sequence, uint32_t start_position, uint32_t end_position) const {
  // Note that the start and end positions should be checked before calling this function
  //canf: kmer_size_ is determined based on the pore model in the Load function above
  int32_t signal_length = end_position - start_position - kmer_size_ + 1;
  assert(signal_length > 0);
  std::vector<float> signal_values;
  signal_values.reserve(signal_length);
  uint32_t mask = ((uint32_t)1 << (2 * kmer_size_)) - 1;
  //canf: These are 2-bit encoded values of k-mers. k-mers are 5/6-mer long so 2^(2*k) many possible values.
  //canf: This only fetches the signal value for the first k-mer. The rest is fetched below in the for loop
  uint32_t hash_value = GenerateSeedFromSequence(sequence, start_position + signal_length, start_position, kmer_size_);
  //canf: Since we have the mean (and stdv) values of each k-mer, we can ask for the signal value of that certain k-mer based on the model.
  //canf: This model should not matter for overlap detection, but it is important for mapping to a reference genome
  signal_values.emplace_back(pore_models_[hash_value].level_mean);
  fprintf(stderr, "%f ", pore_models_[hash_value].level_mean);
//   fprintf(stderr, "%u ", hash_value);
  for (uint32_t position = start_position + 1; position < end_position - kmer_size_ + 1; ++position) {
    uint8_t current_base = CharToUint8(sequence[position + kmer_size_]);
    if (current_base < 4) {
      hash_value = ((hash_value << 2) | current_base) & mask;
    } else {
      hash_value = (hash_value << 2) & mask;
    }

    signal_values.emplace_back(pore_models_[hash_value].level_mean);
	if(position%1000 == 0) fprintf(stderr, "%f ", pore_models_[hash_value].level_mean);
	// if(position < 10000) fprintf(stderr, "%u ", hash_value);
  }

  fprintf(stderr, "\n\n");
  return signal_values;
}
}  // namespace sigmap
