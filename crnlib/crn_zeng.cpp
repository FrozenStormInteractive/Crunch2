// File: crn_zeng.cpp
// See Copyright Notice and license at the end of inc/crnlib.h
// Modified Zeng's technique for codebook/palette reordering
// Evaluation of some reordering techniques for image VQ index compression, António R. C. Paiva , O J. Pinho
// http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.88.7221
#include "crn_core.h"
#include "crn_zeng.h"
#include "crn_sparse_array.h"
#include <deque>

namespace crnlib {
void update_hist(hist_type& hist, int i, int j, int n) {
  if (i == j)
    return;

  if ((i != -1) && (j != -1) && (i < j)) {
    CRNLIB_ASSERT((i >= 0) && (i < (int)n));
    CRNLIB_ASSERT((j >= 0) && (j < (int)n));

    uint index = i * n + j;

#if CRN_ZENG_USE_SPARSE_ARRAY
    uint freq = hist[index];
    freq++;
    hist.set(index, freq);
#else
    hist[index]++;
#endif
  }
}

static inline uint read_hist(hist_type& hist, int i, int j, int n) {
  if (i > j)
    utils::swap(i, j);

  return hist[i * n + j];
}

void create_zeng_reorder_table(uint n, hist_type& xhist, crnlib::vector<uint>& remap_table, zeng_similarity_func pFunc, void* pContext, float similarity_func_weight) {
  CRNLIB_ASSERT(n > 0);
  CRNLIB_ASSERT_CLOSED_RANGE(similarity_func_weight, 0.0f, 1.0f);

  remap_table.clear();
  remap_table.resize(n);

  const uint t = n * n;

  uint max_freq = 0;
  uint max_index = 0;
  for (uint i = 0; i < t; i++) {
    if (xhist[i] > max_freq) {
      max_freq = xhist[i];
      max_index = i;
    }
  }

  uint x = max_index / n;
  uint y = max_index % n;

  crnlib::vector<uint16> values_chosen;
  values_chosen.reserve(n);

  values_chosen.push_back(static_cast<uint16>(x));
  values_chosen.push_back(static_cast<uint16>(y));

  crnlib::vector<uint16> values_remaining;
  if (n > 2)
    values_remaining.reserve(n - 2);
  for (uint i = 0; i < n; i++)
    if ((i != x) && (i != y))
      values_remaining.push_back(static_cast<uint16>(i));

  crnlib::vector<uint> total_freq_to_chosen_values(n);
  for (uint i = 0; i < values_remaining.size(); i++) {
    uint u = values_remaining[i];

    uint total_freq = 0;

    for (uint j = 0; j < values_chosen.size(); j++) {
      uint l = values_chosen[j];

      total_freq += read_hist(xhist, u, l, n);  //[u * n + l];
    }

    total_freq_to_chosen_values[u] = total_freq;
  }

  while (!values_remaining.empty()) {
    double best_freq = 0;
    uint best_i = 0;

    for (uint i = 0; i < values_remaining.size(); i++) {
      uint u = values_remaining[i];
      double total_freq = total_freq_to_chosen_values[u];

      if (pFunc) {
        float weight = math::maximum<float>(
            (*pFunc)(u, values_chosen.front(), pContext),
            (*pFunc)(u, values_chosen.back(), pContext));

        CRNLIB_ASSERT_CLOSED_RANGE(weight, 0.0f, 1.0f);

        weight = math::lerp(1.0f - similarity_func_weight, 1.0f + similarity_func_weight, weight);

        total_freq = (total_freq + 1.0f) * weight;
      }

      if (total_freq > best_freq) {
        best_freq = total_freq;
        best_i = i;
      }
    }

    const uint u = values_remaining[best_i];

    float side = 0;
    int left_freq = 0;
    int right_freq = 0;

    for (uint j = 0; j < values_chosen.size(); j++) {
      const uint l = values_chosen[j];

      int freq = read_hist(xhist, u, l, n);  //[u * n + l];
      int scale = (values_chosen.size() + 1 - 2 * (j + 1));

      side = side + (float)(scale * freq);

      if (scale < 0)
        right_freq += -scale * freq;
      else
        left_freq += scale * freq;
    }

    if (pFunc) {
      float weight_left = (*pFunc)(u, values_chosen.front(), pContext);
      float weight_right = (*pFunc)(u, values_chosen.back(), pContext);

      weight_left = math::lerp(1.0f - similarity_func_weight, 1.0f + similarity_func_weight, weight_left);
      weight_right = math::lerp(1.0f - similarity_func_weight, 1.0f + similarity_func_weight, weight_right);

      side = weight_left * left_freq - weight_right * right_freq;
    }

    if (side > 0)
      values_chosen.push_front(static_cast<uint16>(u));
    else
      values_chosen.push_back(static_cast<uint16>(u));

    values_remaining.erase(values_remaining.begin() + best_i);

    for (uint i = 0; i < values_remaining.size(); i++) {
      const uint r = values_remaining[i];

      total_freq_to_chosen_values[r] += read_hist(xhist, r, u, n);  //[r * n + u];
    }
  }

  for (uint i = 0; i < n; i++) {
    uint v = values_chosen[i];
    remap_table[v] = i;
  }
}

}  // namespace crnlib
