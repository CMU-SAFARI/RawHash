#pragma once

#include <cstdint>
#include <vector>

struct position_pair{
    size_t i; //position in the reference
    size_t j; //position in the read
};

struct alignment_element {
    position_pair position;
    float difference;
};

struct dtw_result {
    float cost;
    std::vector<alignment_element> alignment;
};

float DTW_global(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element = false);
float DTW_global_slow(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element = false);
float DTW_global_diagonalbanded(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, int band_radius, bool exclude_last_element = false);
float DTW_global_slantedbanded(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, int band_radius, bool exclude_last_element = false);
float DTW_global_slantedbanded_antidiagonalwise(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, int band_radius, bool exclude_last_element = false);
float DTW_semiglobal(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element = false);
float DTW_semiglobal_slow(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element = false);
dtw_result DTW_global_tb(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element = false);
dtw_result DTW_semiglobal_tb(const float* a_values, const uint32_t a_length, const float* b_values, const uint32_t b_length, bool exclude_last_element = false);
