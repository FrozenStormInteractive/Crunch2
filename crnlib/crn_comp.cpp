// File: crn_comp.cpp
// See Copyright Notice and license at the end of inc/crnlib.h
#include "crn_core.h"
#include "crn_console.h"
#include "crn_comp.h"
#include "crn_zeng.h"
#include "crn_checksum.h"

#define CRNLIB_CREATE_DEBUG_IMAGES 0
#define CRNLIB_ENABLE_DEBUG_MESSAGES 0

namespace crnlib {
static const uint cEncodingMapNumChunksPerCode = 3;

crn_comp::crn_comp()
    : m_pParams(NULL) {
}

crn_comp::~crn_comp() {
}

float crn_comp::color_endpoint_similarity_func(uint index_a, uint index_b, void* pContext) {
  dxt_hc& hvq = *static_cast<dxt_hc*>(pContext);

  uint endpoint_a = hvq.get_color_endpoint(index_a);
  uint endpoint_b = hvq.get_color_endpoint(index_b);

  color_quad_u8 a[2];
  a[0] = dxt1_block::unpack_color((uint16)(endpoint_a & 0xFFFF), true);
  a[1] = dxt1_block::unpack_color((uint16)((endpoint_a >> 16) & 0xFFFF), true);

  color_quad_u8 b[2];
  b[0] = dxt1_block::unpack_color((uint16)(endpoint_b & 0xFFFF), true);
  b[1] = dxt1_block::unpack_color((uint16)((endpoint_b >> 16) & 0xFFFF), true);

  uint total_error = color::elucidian_distance(a[0], b[0], false) + color::elucidian_distance(a[1], b[1], false);

  float weight = 1.0f - math::clamp(total_error * 1.0f / 8000.0f, 0.0f, 1.0f);
  return weight;
}

float crn_comp::alpha_endpoint_similarity_func(uint index_a, uint index_b, void* pContext) {
  dxt_hc& hvq = *static_cast<dxt_hc*>(pContext);

  uint endpoint_a = hvq.get_alpha_endpoint(index_a);
  int endpoint_a_lo = dxt5_block::unpack_endpoint(endpoint_a, 0);
  int endpoint_a_hi = dxt5_block::unpack_endpoint(endpoint_a, 1);

  uint endpoint_b = hvq.get_alpha_endpoint(index_b);
  int endpoint_b_lo = dxt5_block::unpack_endpoint(endpoint_b, 0);
  int endpoint_b_hi = dxt5_block::unpack_endpoint(endpoint_b, 1);

  int total_error = math::square(endpoint_a_lo - endpoint_b_lo) + math::square(endpoint_a_hi - endpoint_b_hi);

  float weight = 1.0f - math::clamp(total_error * 1.0f / 256.0f, 0.0f, 1.0f);
  return weight;
}

void crn_comp::sort_color_endpoint_codebook(crnlib::vector<uint>& remapping, const crnlib::vector<uint>& endpoints) {
  remapping.resize(endpoints.size());

  uint lowest_energy = UINT_MAX;
  uint lowest_energy_index = 0;

  for (uint i = 0; i < endpoints.size(); i++) {
    color_quad_u8 a(dxt1_block::unpack_color(static_cast<uint16>(endpoints[i] & 0xFFFF), true));
    color_quad_u8 b(dxt1_block::unpack_color(static_cast<uint16>((endpoints[i] >> 16) & 0xFFFF), true));

    uint total = a.r + a.g + a.b + b.r + b.g + b.b;

    if (total < lowest_energy) {
      lowest_energy = total;
      lowest_energy_index = i;
    }
  }

  uint cur_index = lowest_energy_index;

  crnlib::vector<bool> chosen_flags(endpoints.size());

  uint n = 0;
  for (;;) {
    chosen_flags[cur_index] = true;

    remapping[cur_index] = n;
    n++;
    if (n == endpoints.size())
      break;

    uint lowest_error = UINT_MAX;
    uint lowest_error_index = 0;

    color_quad_u8 a(dxt1_block::unpack_endpoint(endpoints[cur_index], 0, true));
    color_quad_u8 b(dxt1_block::unpack_endpoint(endpoints[cur_index], 1, true));

    for (uint i = 0; i < endpoints.size(); i++) {
      if (chosen_flags[i])
        continue;

      color_quad_u8 c(dxt1_block::unpack_endpoint(endpoints[i], 0, true));
      color_quad_u8 d(dxt1_block::unpack_endpoint(endpoints[i], 1, true));

      uint total = color::elucidian_distance(a, c, false) + color::elucidian_distance(b, d, false);

      if (total < lowest_error) {
        lowest_error = total;
        lowest_error_index = i;
      }
    }

    cur_index = lowest_error_index;
  }
}

void crn_comp::sort_alpha_endpoint_codebook(crnlib::vector<uint>& remapping, const crnlib::vector<uint>& endpoints) {
  remapping.resize(endpoints.size());

  uint lowest_energy = UINT_MAX;
  uint lowest_energy_index = 0;

  for (uint i = 0; i < endpoints.size(); i++) {
    uint a = dxt5_block::unpack_endpoint(endpoints[i], 0);
    uint b = dxt5_block::unpack_endpoint(endpoints[i], 1);

    uint total = a + b;

    if (total < lowest_energy) {
      lowest_energy = total;
      lowest_energy_index = i;
    }
  }

  uint cur_index = lowest_energy_index;

  crnlib::vector<bool> chosen_flags(endpoints.size());

  uint n = 0;
  for (;;) {
    chosen_flags[cur_index] = true;

    remapping[cur_index] = n;
    n++;
    if (n == endpoints.size())
      break;

    uint lowest_error = UINT_MAX;
    uint lowest_error_index = 0;

    const int a = dxt5_block::unpack_endpoint(endpoints[cur_index], 0);
    const int b = dxt5_block::unpack_endpoint(endpoints[cur_index], 1);

    for (uint i = 0; i < endpoints.size(); i++) {
      if (chosen_flags[i])
        continue;

      const int c = dxt5_block::unpack_endpoint(endpoints[i], 0);
      const int d = dxt5_block::unpack_endpoint(endpoints[i], 1);

      uint total = math::square(a - c) + math::square(b - d);

      if (total < lowest_error) {
        lowest_error = total;
        lowest_error_index = i;
      }
    }

    cur_index = lowest_error_index;
  }
}

// The indices are only used for statistical purposes.
bool crn_comp::pack_color_endpoints(
    crnlib::vector<uint8>& data,
    const crnlib::vector<uint>& remapping,
    uint trial_index) {
  trial_index;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("pack_color_endpoints: %u", trial_index);
#endif

  crnlib::vector<uint> remapped_endpoints(m_hvq.get_color_endpoint_codebook_size());

  for (uint i = 0; i < m_hvq.get_color_endpoint_codebook_size(); i++)
    remapped_endpoints[remapping[i]] = m_hvq.get_color_endpoint(i);

  const uint component_limits[6] = {31, 63, 31, 31, 63, 31};

  symbol_histogram hist[2];
  hist[0].resize(32);
  hist[1].resize(64);

#if CRNLIB_CREATE_DEBUG_IMAGES
  image_u8 endpoint_image(2, m_hvq.get_color_endpoint_codebook_size());
  image_u8 endpoint_residual_image(2, m_hvq.get_color_endpoint_codebook_size());
#endif

  crnlib::vector<uint> residual_syms;
  residual_syms.reserve(m_hvq.get_color_endpoint_codebook_size() * 2 * 3);

  color_quad_u8 prev[2];
  prev[0].clear();
  prev[1].clear();

  int total_residuals = 0;

  for (uint endpoint_index = 0; endpoint_index < m_hvq.get_color_endpoint_codebook_size(); endpoint_index++) {
    const uint endpoint = remapped_endpoints[endpoint_index];

    color_quad_u8 cur[2];
    cur[0] = dxt1_block::unpack_color((uint16)(endpoint & 0xFFFF), false);
    cur[1] = dxt1_block::unpack_color((uint16)((endpoint >> 16) & 0xFFFF), false);

#if CRNLIB_CREATE_DEBUG_IMAGES
    endpoint_image(0, endpoint_index) = dxt1_block::unpack_color((uint16)(endpoint & 0xFFFF), true);
    endpoint_image(1, endpoint_index) = dxt1_block::unpack_color((uint16)((endpoint >> 16) & 0xFFFF), true);
#endif

    for (uint j = 0; j < 2; j++) {
      for (uint k = 0; k < 3; k++) {
        int delta = cur[j][k] - prev[j][k];
        total_residuals += delta * delta;

        int sym = delta & component_limits[j * 3 + k];
        int table = (k == 1) ? 1 : 0;

        hist[table].inc_freq(sym);

        residual_syms.push_back(sym);

#if CRNLIB_CREATE_DEBUG_IMAGES
        endpoint_residual_image(j, endpoint_index)[k] = static_cast<uint8>(sym);
#endif
      }
    }

    prev[0] = cur[0];
    prev[1] = cur[1];
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("Total endpoint residuals: %i", total_residuals);
#endif

#if CRNLIB_CREATE_DEBUG_IMAGES
  image_utils::write_to_file(dynamic_string(cVarArg, "color_endpoint_residuals_%u.tga", trial_index).get_ptr(), endpoint_residual_image);
  image_utils::write_to_file(dynamic_string(cVarArg, "color_endpoints_%u.tga", trial_index).get_ptr(), endpoint_image);
#endif

  static_huffman_data_model residual_dm[2];

  symbol_codec codec;
  codec.start_encoding(1024 * 1024);

  // Transmit residuals
  for (uint i = 0; i < 2; i++) {
    if (!residual_dm[i].init(true, hist[i], 15))
      return false;

    if (!codec.encode_transmit_static_huffman_data_model(residual_dm[i], false))
      return false;
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("Wrote %u bits for color endpoint residual Huffman tables", codec.encode_get_total_bits_written());
#endif

  uint start_bits = codec.encode_get_total_bits_written();
  start_bits;

  for (uint i = 0; i < residual_syms.size(); i++) {
    const uint sym = residual_syms[i];
    const uint table = ((i % 3) == 1) ? 1 : 0;
    codec.encode(sym, residual_dm[table]);
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("Wrote %u bits for color endpoint residuals", codec.encode_get_total_bits_written() - start_bits);
#endif

  codec.stop_encoding(false);

  data.swap(codec.get_encoding_buf());

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging) {
    console::debug("Wrote a total of %u bits for color endpoint codebook", codec.encode_get_total_bits_written());

    console::debug("Wrote %f bits per each color endpoint", data.size() * 8.0f / m_hvq.get_color_endpoint_codebook_size());
  }
#endif

  return true;
}

// The indices are only used for statistical purposes.
bool crn_comp::pack_alpha_endpoints(
    crnlib::vector<uint8>& data,
    const crnlib::vector<uint>& remapping,
    uint trial_index) {
  trial_index;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("pack_alpha_endpoints: %u", trial_index);
#endif

  crnlib::vector<uint> remapped_endpoints(m_hvq.get_alpha_endpoint_codebook_size());

  for (uint i = 0; i < m_hvq.get_alpha_endpoint_codebook_size(); i++)
    remapped_endpoints[remapping[i]] = m_hvq.get_alpha_endpoint(i);

  symbol_histogram hist;
  hist.resize(256);

#if CRNLIB_CREATE_DEBUG_IMAGES
  image_u8 endpoint_image(2, m_hvq.get_alpha_endpoint_codebook_size());
  image_u8 endpoint_residual_image(2, m_hvq.get_alpha_endpoint_codebook_size());
#endif

  crnlib::vector<uint> residual_syms;
  residual_syms.reserve(m_hvq.get_alpha_endpoint_codebook_size() * 2 * 3);

  uint prev[2];
  utils::zero_object(prev);

  int total_residuals = 0;

  for (uint endpoint_index = 0; endpoint_index < m_hvq.get_alpha_endpoint_codebook_size(); endpoint_index++) {
    const uint endpoint = remapped_endpoints[endpoint_index];

    uint cur[2];
    cur[0] = dxt5_block::unpack_endpoint(endpoint, 0);
    cur[1] = dxt5_block::unpack_endpoint(endpoint, 1);

#if CRNLIB_CREATE_DEBUG_IMAGES
    endpoint_image(0, endpoint_index) = cur[0];
    endpoint_image(1, endpoint_index) = cur[1];
#endif

    for (uint j = 0; j < 2; j++) {
      int delta = cur[j] - prev[j];
      total_residuals += delta * delta;

      int sym = delta & 255;

      hist.inc_freq(sym);

      residual_syms.push_back(sym);

#if CRNLIB_CREATE_DEBUG_IMAGES
      endpoint_residual_image(j, endpoint_index) = static_cast<uint8>(sym);
#endif
    }

    prev[0] = cur[0];
    prev[1] = cur[1];
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("Total endpoint residuals: %i", total_residuals);
#endif

#if CRNLIB_CREATE_DEBUG_IMAGES
  image_utils::write_to_file(dynamic_string(cVarArg, "alpha_endpoint_residuals_%u.tga", trial_index).get_ptr(), endpoint_residual_image);
  image_utils::write_to_file(dynamic_string(cVarArg, "alpha_endpoints_%u.tga", trial_index).get_ptr(), endpoint_image);
#endif

  static_huffman_data_model residual_dm;

  symbol_codec codec;
  codec.start_encoding(1024 * 1024);

  // Transmit residuals
  if (!residual_dm.init(true, hist, 15))
    return false;

  if (!codec.encode_transmit_static_huffman_data_model(residual_dm, false))
    return false;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("Wrote %u bits for alpha endpoint residual Huffman tables", codec.encode_get_total_bits_written());
#endif

  uint start_bits = codec.encode_get_total_bits_written();
  start_bits;

  for (uint i = 0; i < residual_syms.size(); i++) {
    const uint sym = residual_syms[i];
    codec.encode(sym, residual_dm);
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("Wrote %u bits for alpha endpoint residuals", codec.encode_get_total_bits_written() - start_bits);
#endif

  codec.stop_encoding(false);

  data.swap(codec.get_encoding_buf());

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging) {
    console::debug("Wrote a total of %u bits for alpha endpoint codebook", codec.encode_get_total_bits_written());

    console::debug("Wrote %f bits per each alpha endpoint", data.size() * 8.0f / m_hvq.get_alpha_endpoint_codebook_size());
  }
#endif

  return true;
}

float crn_comp::color_selector_similarity_func(uint index_a, uint index_b, void* pContext) {
  const crnlib::vector<dxt_hc::selectors>& selectors = *static_cast<const crnlib::vector<dxt_hc::selectors>*>(pContext);

  const dxt_hc::selectors& selectors_a = selectors[index_a];
  const dxt_hc::selectors& selectors_b = selectors[index_b];

  int total = 0;
  for (uint i = 0; i < 16; i++) {
    int a = g_dxt1_to_linear[selectors_a.get_by_index(i)];
    int b = g_dxt1_to_linear[selectors_b.get_by_index(i)];

    int delta = a - b;
    total += delta * delta;
  }

  float weight = 1.0f - math::clamp(total * 1.0f / 20.0f, 0.0f, 1.0f);
  return weight;
}

float crn_comp::alpha_selector_similarity_func(uint index_a, uint index_b, void* pContext) {
  const crnlib::vector<dxt_hc::selectors>& selectors = *static_cast<const crnlib::vector<dxt_hc::selectors>*>(pContext);

  const dxt_hc::selectors& selectors_a = selectors[index_a];
  const dxt_hc::selectors& selectors_b = selectors[index_b];

  int total = 0;
  for (uint i = 0; i < 16; i++) {
    int a = g_dxt5_to_linear[selectors_a.get_by_index(i)];
    int b = g_dxt5_to_linear[selectors_b.get_by_index(i)];

    int delta = a - b;
    total += delta * delta;
  }

  float weight = 1.0f - math::clamp(total * 1.0f / 100.0f, 0.0f, 1.0f);
  return weight;
}

void crn_comp::sort_selector_codebook(crnlib::vector<uint>& remapping, const crnlib::vector<dxt_hc::selectors>& selectors, const uint8* pTo_linear) {
  remapping.resize(selectors.size());

  uint lowest_energy = UINT_MAX;
  uint lowest_energy_index = 0;

  for (uint i = 0; i < selectors.size(); i++) {
    uint total = 0;
    for (uint j = 0; j < 16; j++) {
      int a = pTo_linear[selectors[i].get_by_index(j)];

      total += a * a;
    }

    if (total < lowest_energy) {
      lowest_energy = total;
      lowest_energy_index = i;
    }
  }

  uint cur_index = lowest_energy_index;

  crnlib::vector<bool> chosen_flags(selectors.size());

  uint n = 0;
  for (;;) {
    chosen_flags[cur_index] = true;

    remapping[cur_index] = n;
    n++;
    if (n == selectors.size())
      break;

    uint lowest_error = UINT_MAX;
    uint lowest_error_index = 0;

    for (uint i = 0; i < selectors.size(); i++) {
      if (chosen_flags[i])
        continue;

      uint total = 0;
      for (uint j = 0; j < 16; j++) {
        int a = pTo_linear[selectors[cur_index].get_by_index(j)];
        int b = pTo_linear[selectors[i].get_by_index(j)];

        int delta = a - b;
        total += delta * delta;
      }

      if (total < lowest_error) {
        lowest_error = total;
        lowest_error_index = i;
      }
    }

    cur_index = lowest_error_index;
  }
}

// The indices are only used for statistical purposes.
bool crn_comp::pack_selectors(
    crnlib::vector<uint8>& packed_data,
    const crnlib::vector<dxt_hc::selectors>& selectors,
    const crnlib::vector<uint>& remapping,
    uint max_selector_value,
    const uint8* pTo_linear,
    uint trial_index) {
  trial_index;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("pack_selectors: %u", trial_index);
#endif

  crnlib::vector<dxt_hc::selectors> remapped_selectors(selectors.size());

  for (uint i = 0; i < selectors.size(); i++)
    remapped_selectors[remapping[i]] = selectors[i];

#if CRNLIB_CREATE_DEBUG_IMAGES
  image_u8 residual_image(16, selectors.size());
  ;
  image_u8 selector_image(16, selectors.size());
  ;
#endif

  crnlib::vector<uint> residual_syms;
  residual_syms.reserve(selectors.size() * 8);

  const uint num_baised_selector_values = (max_selector_value * 2 + 1);
  symbol_histogram hist(num_baised_selector_values * num_baised_selector_values);

  dxt_hc::selectors prev_selectors;
  utils::zero_object(prev_selectors);
  int total_residuals = 0;
  for (uint selector_index = 0; selector_index < selectors.size(); selector_index++) {
    const dxt_hc::selectors& s = remapped_selectors[selector_index];

    uint prev_sym = 0;
    for (uint i = 0; i < 16; i++) {
      int p = pTo_linear[crnlib_assert_range_incl<uint>(prev_selectors.get_by_index(i), max_selector_value)];

      int r = pTo_linear[crnlib_assert_range_incl<uint>(s.get_by_index(i), max_selector_value)] - p;

      total_residuals += r * r;

      uint sym = r + max_selector_value;

      CRNLIB_ASSERT(sym < num_baised_selector_values);
      if (i & 1) {
        uint paired_sym = (sym * num_baised_selector_values) + prev_sym;
        residual_syms.push_back(paired_sym);
        hist.inc_freq(paired_sym);
      } else
        prev_sym = sym;

#if CRNLIB_CREATE_DEBUG_IMAGES
      selector_image(i, selector_index) = (pTo_linear[crnlib_assert_range_incl<uint>(s.get_by_index(i), max_selector_value)] * 255) / max_selector_value;
      residual_image(i, selector_index) = sym;
#endif
    }

    prev_selectors = s;
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("Total selector endpoint residuals: %u", total_residuals);
#endif

#if CRNLIB_CREATE_DEBUG_IMAGES
  image_utils::write_to_file(dynamic_string(cVarArg, "selectors_%u_%u.tga", trial_index, max_selector_value).get_ptr(), selector_image);
  image_utils::write_to_file(dynamic_string(cVarArg, "selector_residuals_%u_%u.tga", trial_index, max_selector_value).get_ptr(), residual_image);
#endif

  static_huffman_data_model residual_dm;

  symbol_codec codec;
  codec.start_encoding(1024 * 1024);

  // Transmit residuals
  if (!residual_dm.init(true, hist, 15))
    return false;

  if (!codec.encode_transmit_static_huffman_data_model(residual_dm, false))
    return false;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("Wrote %u bits for selector residual Huffman tables", codec.encode_get_total_bits_written());
#endif

  uint start_bits = codec.encode_get_total_bits_written();
  start_bits;

  for (uint i = 0; i < residual_syms.size(); i++) {
    const uint sym = residual_syms[i];
    codec.encode(sym, residual_dm);
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("Wrote %u bits for selector residuals", codec.encode_get_total_bits_written() - start_bits);
#endif

  codec.stop_encoding(false);

  packed_data.swap(codec.get_encoding_buf());

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging) {
    console::debug("Wrote a total of %u bits for selector codebook", codec.encode_get_total_bits_written());

    console::debug("Wrote %f bits per each selector codebook entry", packed_data.size() * 8.0f / selectors.size());
  }
#endif

  return true;
}

bool crn_comp::pack_chunks(
    uint group,
    bool clear_histograms,
    symbol_codec* pCodec,
    const crnlib::vector<uint>* pColor_endpoint_remap,
    const crnlib::vector<uint>* pColor_selector_remap,
    const crnlib::vector<uint>* pAlpha_endpoint_remap,
    const crnlib::vector<uint>* pAlpha_selector_remap) {
  uint first_chunk = m_mip_groups[group].m_first_chunk;
  uint num_chunks = m_mip_groups[group].m_num_chunks;
  uint chunk_width = m_mip_groups[group].m_chunk_width;

  if (!pCodec) {
    m_chunk_encoding_hist.resize(256);
    if (clear_histograms)
      m_chunk_encoding_hist.set_all(0);

    if (pColor_endpoint_remap) {
      CRNLIB_ASSERT(pColor_endpoint_remap->size() == m_hvq.get_color_endpoint_codebook_size());
      m_endpoint_index_hist[0].resize(pColor_endpoint_remap->size());
      if (clear_histograms)
        m_endpoint_index_hist[0].set_all(0);
    }

    if (pColor_selector_remap) {
      CRNLIB_ASSERT(pColor_selector_remap->size() == m_hvq.get_color_selector_codebook_size());
      m_selector_index_hist[0].resize(pColor_selector_remap->size());
      if (clear_histograms)
        m_selector_index_hist[0].set_all(0);
    }

    if (pAlpha_endpoint_remap) {
      CRNLIB_ASSERT(pAlpha_endpoint_remap->size() == m_hvq.get_alpha_endpoint_codebook_size());
      m_endpoint_index_hist[1].resize(pAlpha_endpoint_remap->size());
      if (clear_histograms)
        m_endpoint_index_hist[1].set_all(0);
    }

    if (pAlpha_selector_remap) {
      CRNLIB_ASSERT(pAlpha_selector_remap->size() == m_hvq.get_alpha_selector_codebook_size());
      m_selector_index_hist[1].resize(pAlpha_selector_remap->size());
      if (clear_histograms)
        m_selector_index_hist[1].set_all(0);
    }
  }

  uint endpoint_index[cNumComps] = {};
  const crnlib::vector<uint>* endpoint_remap[cNumComps] = {};
  const crnlib::vector<uint>* selector_remap[cNumComps] = {};
  for (uint c = 0; c < cNumComps; c++) {
    if (m_has_comp[c]) {
      endpoint_remap[c] = c ? pAlpha_endpoint_remap : pColor_endpoint_remap;
      selector_remap[c] = c ? pAlpha_selector_remap : pColor_selector_remap;
    }
  }

  for (uint by = 0, block_width = chunk_width << 1, b = first_chunk << 2, bEnd = b + (num_chunks << 2); b < bEnd; by++) {
    for (uint bx = 0; bx < block_width; bx++, b++) {
      if (!(by & 1) && !(bx & 1)) {
        uint8 reference_group = m_endpoint_indices[b].reference | m_endpoint_indices[b + block_width].reference << 2 |
          m_endpoint_indices[b + 1].reference << 4 | m_endpoint_indices[b + block_width + 1].reference << 6;
        if (pCodec)
          pCodec->encode(reference_group, m_reference_encoding_dm);
        else
          m_chunk_encoding_hist.inc_freq(reference_group);
      }
      for (uint c = 0; c < cNumComps; c++) {
        if (endpoint_remap[c]) {
          uint index = (*endpoint_remap[c])[m_endpoint_indices[b].component[c]];
          if (!m_endpoint_indices[b].reference) {
            int sym = index - endpoint_index[c];
            if (sym < 0)
              sym += endpoint_remap[c]->size();
            if (!pCodec)
              m_endpoint_index_hist[c ? 1 : 0].inc_freq(sym);
            else
              pCodec->encode(sym, m_endpoint_index_dm[c ? 1 : 0]);
          }
          endpoint_index[c] = index;
        }
      }
      for (uint c = 0; c < cNumComps; c++) {
        if (selector_remap[c]) {
          uint index = (*selector_remap[c])[m_selector_indices[b].component[c]];
          if (!pCodec)
            m_selector_index_hist[c ? 1 : 0].inc_freq(index);
          else
            pCodec->encode(index, m_selector_index_dm[c ? 1 : 0]);
        }
      }
    }
  }
  return true;
}

void crn_comp::append_vec(crnlib::vector<uint8>& a, const void* p, uint size) {
  if (size) {
    uint ofs = a.size();
    a.resize(ofs + size);

    memcpy(&a[ofs], p, size);
  }
}

void crn_comp::append_vec(crnlib::vector<uint8>& a, const crnlib::vector<uint8>& b) {
  if (!b.empty()) {
    uint ofs = a.size();
    a.resize(ofs + b.size());

    memcpy(&a[ofs], &b[0], b.size());
  }
}

#if 0
   bool crn_comp::init_chunk_encoding_dm()
   {
      symbol_histogram hist(1 << (3 * cEncodingMapNumChunksPerCode));

      for (uint chunk_index = 0; chunk_index < m_hvq.get_num_chunks(); chunk_index += cEncodingMapNumChunksPerCode)
      {
         uint index = 0;
         for (uint i = 0; i < cEncodingMapNumChunksPerCode; i++)
         {
            if ((chunk_index + i) >= m_hvq.get_num_chunks())
               break;
            const dxt_hc::chunk_encoding& encoding = m_hvq.get_chunk_encoding(chunk_index + i);

            index |= (encoding.m_encoding_index << (i * 3));
         }

         hist.inc_freq(index);
      }

      if (!m_chunk_encoding_dm.init(true, hist, 16))
         return false;

      return true;
   }
#endif

bool crn_comp::alias_images() {
  for (uint face_index = 0; face_index < m_pParams->m_faces; face_index++) {
    for (uint level_index = 0; level_index < m_pParams->m_levels; level_index++) {
      const uint width = math::maximum(1U, m_pParams->m_width >> level_index);
      const uint height = math::maximum(1U, m_pParams->m_height >> level_index);

      if (!m_pParams->m_pImages[face_index][level_index])
        return false;

      m_images[face_index][level_index].alias((color_quad_u8*)m_pParams->m_pImages[face_index][level_index], width, height);
    }
  }

  image_utils::conversion_type conv_type = image_utils::get_image_conversion_type_from_crn_format((crn_format)m_pParams->m_format);
  if (conv_type != image_utils::cConversion_Invalid) {
    for (uint face_index = 0; face_index < m_pParams->m_faces; face_index++) {
      for (uint level_index = 0; level_index < m_pParams->m_levels; level_index++) {
        image_u8 cooked_image(m_images[face_index][level_index]);

        image_utils::convert_image(cooked_image, conv_type);

        m_images[face_index][level_index].swap(cooked_image);
      }
    }
  }

  m_mip_groups.clear();
  m_mip_groups.resize(m_pParams->m_levels);

  utils::zero_object(m_levels);

  uint mip_group = 0;
  uint chunk_index = 0;
  uint mip_group_chunk_index = 0;
  (void)mip_group_chunk_index;
  for (uint level_index = 0; level_index < m_pParams->m_levels; level_index++) {
    const uint width = math::maximum(1U, m_pParams->m_width >> level_index);
    const uint height = math::maximum(1U, m_pParams->m_height >> level_index);
    const uint chunk_width = math::align_up_value(width, cChunkPixelWidth) / cChunkPixelWidth;
    const uint chunk_height = math::align_up_value(height, cChunkPixelHeight) / cChunkPixelHeight;
    const uint num_chunks = m_pParams->m_faces * chunk_width * chunk_height;

    m_mip_groups[mip_group].m_first_chunk = chunk_index;
    mip_group_chunk_index = 0;

    m_mip_groups[mip_group].m_num_chunks += num_chunks;
    m_mip_groups[mip_group].m_chunk_width = chunk_width;

    m_levels[level_index].m_width = width;
    m_levels[level_index].m_height = height;
    m_levels[level_index].m_chunk_width = chunk_width;
    m_levels[level_index].m_chunk_height = chunk_height;
    m_levels[level_index].m_first_chunk = chunk_index;
    m_levels[level_index].m_num_chunks = num_chunks;
    m_levels[level_index].m_group_index = mip_group;
    m_levels[level_index].m_group_first_chunk = 0;

    chunk_index += num_chunks;

    mip_group++;
  }

  m_total_chunks = chunk_index;

  return true;
}

void crn_comp::clear() {
  m_pParams = NULL;

  for (uint f = 0; f < cCRNMaxFaces; f++)
    for (uint l = 0; l < cCRNMaxLevels; l++)
      m_images[f][l].clear();

  utils::zero_object(m_levels);

  m_mip_groups.clear();

  utils::zero_object(m_has_comp);

  m_endpoint_indices.clear();
  m_selector_indices.clear();

  m_total_chunks = 0;

  utils::zero_object(m_crn_header);

  m_comp_data.clear();

  m_hvq.clear();

  m_chunk_encoding_hist.clear();
  m_reference_encoding_dm.clear();
  for (uint i = 0; i < 2; i++) {
    m_endpoint_index_hist[i].clear();
    m_endpoint_index_dm[i].clear();
    m_selector_index_hist[i].clear();
    m_selector_index_dm[i].clear();
  }

  for (uint i = 0; i < cCRNMaxLevels; i++)
    m_packed_chunks[i].clear();

  m_packed_data_models.clear();

  m_packed_color_endpoints.clear();
  m_packed_color_selectors.clear();
  m_packed_alpha_endpoints.clear();
  m_packed_alpha_selectors.clear();
}

bool crn_comp::quantize_images() {
  dxt_hc::params params;

  params.m_adaptive_tile_alpha_psnr_derating = m_pParams->m_crn_adaptive_tile_alpha_psnr_derating;
  params.m_adaptive_tile_color_psnr_derating = m_pParams->m_crn_adaptive_tile_color_psnr_derating;

  if (m_pParams->m_flags & cCRNCompFlagManualPaletteSizes) {
    params.m_color_endpoint_codebook_size = math::clamp<int>(m_pParams->m_crn_color_endpoint_palette_size, cCRNMinPaletteSize, cCRNMaxPaletteSize);
    params.m_color_selector_codebook_size = math::clamp<int>(m_pParams->m_crn_color_selector_palette_size, cCRNMinPaletteSize, cCRNMaxPaletteSize);
    params.m_alpha_endpoint_codebook_size = math::clamp<int>(m_pParams->m_crn_alpha_endpoint_palette_size, cCRNMinPaletteSize, cCRNMaxPaletteSize);
    params.m_alpha_selector_codebook_size = math::clamp<int>(m_pParams->m_crn_alpha_selector_palette_size, cCRNMinPaletteSize, cCRNMaxPaletteSize);
  } else {
    uint max_codebook_entries = ((m_pParams->m_width + 3) / 4) * ((m_pParams->m_height + 3) / 4);

    max_codebook_entries = math::clamp<uint>(max_codebook_entries, cCRNMinPaletteSize, cCRNMaxPaletteSize);

    float quality = math::clamp<float>((float)m_pParams->m_quality_level / cCRNMaxQualityLevel, 0.0f, 1.0f);
    float color_quality_power_mul = 1.0f;
    float alpha_quality_power_mul = 1.0f;
    if (m_pParams->m_format == cCRNFmtDXT5_CCxY) {
      color_quality_power_mul = 3.5f;
      alpha_quality_power_mul = .35f;
      params.m_adaptive_tile_color_psnr_derating = 5.0f;
    } else if (m_pParams->m_format == cCRNFmtDXT5)
      color_quality_power_mul = .75f;

    float color_endpoint_quality = powf(quality, 1.8f * color_quality_power_mul);
    float color_selector_quality = powf(quality, 1.65f * color_quality_power_mul);
    params.m_color_endpoint_codebook_size = math::clamp<uint>(math::float_to_uint(.5f + math::lerp<float>(math::maximum<float>(64, cCRNMinPaletteSize), (float)max_codebook_entries, color_endpoint_quality)), cCRNMinPaletteSize, cCRNMaxPaletteSize);
    params.m_color_selector_codebook_size = math::clamp<uint>(math::float_to_uint(.5f + math::lerp<float>(math::maximum<float>(96, cCRNMinPaletteSize), (float)max_codebook_entries, color_selector_quality)), cCRNMinPaletteSize, cCRNMaxPaletteSize);

    float alpha_endpoint_quality = powf(quality, 2.1f * alpha_quality_power_mul);
    float alpha_selector_quality = powf(quality, 1.65f * alpha_quality_power_mul);
    params.m_alpha_endpoint_codebook_size = math::clamp<uint>(math::float_to_uint(.5f + math::lerp<float>(math::maximum<float>(24, cCRNMinPaletteSize), (float)max_codebook_entries, alpha_endpoint_quality)), cCRNMinPaletteSize, cCRNMaxPaletteSize);    
    params.m_alpha_selector_codebook_size = math::clamp<uint>(math::float_to_uint(.5f + math::lerp<float>(math::maximum<float>(48, cCRNMinPaletteSize), (float)max_codebook_entries, alpha_selector_quality)), cCRNMinPaletteSize, cCRNMaxPaletteSize);
  }

  if (m_pParams->m_flags & cCRNCompFlagDebugging) {
    console::debug("Color endpoints: %u", params.m_color_endpoint_codebook_size);
    console::debug("Color selectors: %u", params.m_color_selector_codebook_size);
    console::debug("Alpha endpoints: %u", params.m_alpha_endpoint_codebook_size);
    console::debug("Alpha selectors: %u", params.m_alpha_selector_codebook_size);
  }

  params.m_hierarchical = (m_pParams->m_flags & cCRNCompFlagHierarchical) != 0;
  params.m_perceptual = (m_pParams->m_flags & cCRNCompFlagPerceptual) != 0;

  params.m_pProgress_func = m_pParams->m_pProgress_func;
  params.m_pProgress_func_data = m_pParams->m_pProgress_func_data;

  switch (m_pParams->m_format) {
    case cCRNFmtDXT1: {
      params.m_format = cDXT1;
      m_has_comp[cColor] = true;
      break;
    }
    case cCRNFmtDXT3: {
      m_has_comp[cAlpha0] = true;
      return false;
    }
    case cCRNFmtDXT5: {
      params.m_format = cDXT5;
      params.m_alpha_component_indices[0] = m_pParams->m_alpha_component;
      m_has_comp[cColor] = true;
      m_has_comp[cAlpha0] = true;
      break;
    }
    case cCRNFmtDXT5_CCxY: {
      params.m_format = cDXT5;
      params.m_alpha_component_indices[0] = 3;
      m_has_comp[cColor] = true;
      m_has_comp[cAlpha0] = true;
      params.m_perceptual = false;

      //params.m_adaptive_tile_color_alpha_weighting_ratio = 1.0f;
      params.m_adaptive_tile_color_alpha_weighting_ratio = 1.5f;
      break;
    }
    case cCRNFmtDXT5_xGBR:
    case cCRNFmtDXT5_AGBR:
    case cCRNFmtDXT5_xGxR: {
      params.m_format = cDXT5;
      params.m_alpha_component_indices[0] = 3;
      m_has_comp[cColor] = true;
      m_has_comp[cAlpha0] = true;
      params.m_perceptual = false;
      break;
    }
    case cCRNFmtDXN_XY: {
      params.m_format = cDXN_XY;
      params.m_alpha_component_indices[0] = 0;
      params.m_alpha_component_indices[1] = 1;
      m_has_comp[cAlpha0] = true;
      m_has_comp[cAlpha1] = true;
      params.m_perceptual = false;
      break;
    }
    case cCRNFmtDXN_YX: {
      params.m_format = cDXN_YX;
      params.m_alpha_component_indices[0] = 1;
      params.m_alpha_component_indices[1] = 0;
      m_has_comp[cAlpha0] = true;
      m_has_comp[cAlpha1] = true;
      params.m_perceptual = false;
      break;
    }
    case cCRNFmtDXT5A: {
      params.m_format = cDXT5A;
      params.m_alpha_component_indices[0] = m_pParams->m_alpha_component;
      m_has_comp[cAlpha0] = true;
      params.m_perceptual = false;
      break;
    }
    case cCRNFmtETC1: {
      console::warning("crn_comp::quantize_chunks: This class does not support ETC1");
      return false;
    }
    default: {
      return false;
    }
  }
  params.m_debugging = (m_pParams->m_flags & cCRNCompFlagDebugging) != 0;

  params.m_num_levels = m_pParams->m_levels;
  for (uint i = 0; i < m_pParams->m_levels; i++) {
    params.m_levels[i].m_first_block = m_levels[i].m_first_chunk << 2;
    params.m_levels[i].m_num_blocks = m_levels[i].m_num_chunks << 2;
    params.m_levels[i].m_block_width = m_levels[i].m_chunk_width << 1;
    params.m_levels[i].m_weight = math::minimum(12.0f, powf(1.3f, (float)i));
  }
  params.m_num_faces = m_pParams->m_faces;

  params.m_endpoint_indices = &m_endpoint_indices;
  params.m_selector_indices = &m_selector_indices;

  params.m_num_blocks = m_total_chunks << 2;
  params.m_blocks = (color_quad_u8(*)[16])crnlib_malloc(params.m_num_blocks * 16 * sizeof(color_quad_u8));
  for (uint b = 0, level = 0; level < m_pParams->m_levels; level++) {
    for (uint face = 0; face < m_pParams->m_faces; face++) {
      image_u8& image = m_images[face][level];
      uint width = image.get_width();
      uint height = image.get_height();
      uint blockWidth = (width + 7 & ~7) >> 2;
      uint blockHeight = (height + 7 & ~7) >> 2;
      for (uint by = 0; by < blockHeight; by++) {
        for (uint y0 = by << 2, bx = 0; bx < blockWidth; bx++, b++) {
          for (uint t = 0, x0 = bx << 2, dy = 0; dy < 4; dy++) {
            for (uint y = math::minimum<uint>(y0 + dy, height - 1), dx = 0; dx < 4; dx++, t++)
              params.m_blocks[b][t] = image(math::minimum<uint>(x0 + dx, width - 1), y);
          }
        }
      }
    }
  }
  bool result = m_hvq.compress(params, m_task_pool);
  crnlib_free(params.m_blocks);

  return result;
}

struct optimize_color_endpoint_codebook_params {
  struct trial {
    crnlib::vector<uint> remapping;
    crnlib::vector<uint8> packed_endpoints;
    uint total_bits;
  } *m_trial;
  hist_type* m_xhist;
  uint m_iter_index;
  uint m_max_iter_index;
};

void crn_comp::optimize_color_endpoint_codebook_task(uint64 data, void* pData_ptr) {
  data;
  optimize_color_endpoint_codebook_params* pParams = reinterpret_cast<optimize_color_endpoint_codebook_params*>(pData_ptr);
  crnlib::vector<uint>& remapping = pParams->m_trial->remapping;

  if (pParams->m_iter_index == pParams->m_max_iter_index) {
    sort_color_endpoint_codebook(remapping, m_hvq.get_color_endpoint_vec());
  } else {
    create_zeng_reorder_table(
        m_hvq.get_color_endpoint_codebook_size(),
        *pParams->m_xhist,
        remapping,
        pParams->m_iter_index ? color_endpoint_similarity_func : NULL,
        &m_hvq,
        pParams->m_iter_index / static_cast<float>(pParams->m_max_iter_index - 1));
  }

  pack_color_endpoints(pParams->m_trial->packed_endpoints, remapping, pParams->m_iter_index);
  uint total_bits = pParams->m_trial->packed_endpoints.size() << 3;
  uint codebook_size = remapping.size();

  crnlib::vector<uint> hist(codebook_size);
  for (uint group = 0; group < m_mip_groups.size(); group++) {
    for (uint endpoint_index = 0, b = m_mip_groups[group].m_first_chunk << 2, bEnd = b + (m_mip_groups[group].m_num_chunks << 2); b < bEnd; b++) {
      uint index = remapping[m_endpoint_indices[b].component[cColor]];
      if (!m_endpoint_indices[b].reference) {
        int sym = index - endpoint_index;
        hist[sym < 0 ? sym + codebook_size : sym]++;
      }
      endpoint_index = index;
    }
  }

  static_huffman_data_model dm;
  dm.init(true, codebook_size, hist.get_ptr(), 16);
  const uint8* code_sizes = dm.get_code_sizes();
  for (uint s = 0; s < codebook_size; s++)
    total_bits += hist[s] * code_sizes[s];

  symbol_codec codec;
  codec.start_encoding(64 * 1024);
  codec.encode_enable_simulation(true);
  codec.encode_transmit_static_huffman_data_model(dm, false);
  codec.stop_encoding(false);
  total_bits += codec.encode_get_total_bits_written();

  pParams->m_trial->total_bits = total_bits;

  crnlib_delete(pParams);
}

bool crn_comp::optimize_color_endpoint_codebook(crnlib::vector<uint>& remapping) {
  if (m_pParams->m_flags & cCRNCompFlagQuick) {
    remapping.resize(m_hvq.get_color_endpoint_vec().size());
    for (uint i = 0; i < m_hvq.get_color_endpoint_vec().size(); i++)
      remapping[i] = i;

    if (!pack_color_endpoints(m_packed_color_endpoints, remapping, 0))
      return false;

    return true;
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("----- Begin optimization of color endpoint codebook");
#endif

  const uint cMaxEndpointRemapIters = 3;
  optimize_color_endpoint_codebook_params::trial remapping_trial[cMaxEndpointRemapIters + 1];

  uint n = m_hvq.get_color_endpoint_codebook_size();
  hist_type xhist(n * n);
  for (uint b = 1; b < m_endpoint_indices.size(); b++) {
    if (!m_endpoint_indices[b].reference) {
      update_hist(xhist, m_endpoint_indices[b - 1].color, m_endpoint_indices[b].color, n);
      update_hist(xhist, m_endpoint_indices[b].color, m_endpoint_indices[b - 1].color, n);
    }
  }

  for (uint i = 0; i <= cMaxEndpointRemapIters; i++) {
    optimize_color_endpoint_codebook_params* pParams = crnlib_new<optimize_color_endpoint_codebook_params>();
    pParams->m_iter_index = i;
    pParams->m_max_iter_index = cMaxEndpointRemapIters;
    pParams->m_trial = remapping_trial + i;
    pParams->m_xhist = &xhist;
    m_task_pool.queue_object_task(this, &crn_comp::optimize_color_endpoint_codebook_task, 0, pParams);
  }
  m_task_pool.join();

  for (uint best_bits = UINT_MAX, i = 0; i <= cMaxEndpointRemapIters; i++) {
    if (remapping_trial[i].total_bits < best_bits) {
      m_packed_color_endpoints.swap(remapping_trial[i].packed_endpoints);
      remapping.swap(remapping_trial[i].remapping);
      best_bits = remapping_trial[i].total_bits;
    }
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("End optimization of color endpoint codebook");
#endif

  return true;
}

bool crn_comp::optimize_color_selector_codebook(crnlib::vector<uint>& remapping) {
  if (m_pParams->m_flags & cCRNCompFlagQuick) {
    remapping.resize(m_hvq.get_color_selectors_vec().size());
    for (uint i = 0; i < m_hvq.get_color_selectors_vec().size(); i++)
      remapping[i] = i;
  } else {
    sort_selector_codebook(remapping, m_hvq.get_color_selectors_vec(), g_dxt1_to_linear);
  }
  return pack_selectors(m_packed_color_selectors, m_hvq.get_color_selectors_vec(), remapping, 3, g_dxt1_to_linear, 0);
}

struct optimize_alpha_endpoint_codebook_params {
  struct trial {
    crnlib::vector<uint> remapping;
    crnlib::vector<uint8> packed_endpoints;
    uint total_bits;
  } *m_trial;
  hist_type* m_xhist;
  uint m_iter_index;
  uint m_max_iter_index;
};

void crn_comp::optimize_alpha_endpoint_codebook_task(uint64 data, void* pData_ptr) {
  data;
  optimize_alpha_endpoint_codebook_params* pParams = reinterpret_cast<optimize_alpha_endpoint_codebook_params*>(pData_ptr);
  crnlib::vector<uint>& remapping = pParams->m_trial->remapping;

  if (pParams->m_iter_index == pParams->m_max_iter_index) {
    sort_alpha_endpoint_codebook(remapping, m_hvq.get_alpha_endpoint_vec());
  } else {
    create_zeng_reorder_table(
        m_hvq.get_alpha_endpoint_codebook_size(),
        *pParams->m_xhist,
        remapping,
        pParams->m_iter_index ? alpha_endpoint_similarity_func : NULL,
        &m_hvq,
        pParams->m_iter_index / static_cast<float>(pParams->m_max_iter_index - 1));
  }

  pack_alpha_endpoints(pParams->m_trial->packed_endpoints, remapping, pParams->m_iter_index);
  uint total_bits = pParams->m_trial->packed_endpoints.size() << 3;
  uint codebook_size = remapping.size();

  crnlib::vector<uint> hist(codebook_size);
  bool hasAlpha0 = m_has_comp[cAlpha0], hasAlpha1 = m_has_comp[cAlpha1];
  for (uint group = 0; group < m_mip_groups.size(); group++) {
    for (uint index0 = 0, index1 = 0, b = m_mip_groups[group].m_first_chunk << 2, bEnd = b + (m_mip_groups[group].m_num_chunks << 2); b < bEnd; b++) {
      if (hasAlpha0) {
        uint index = remapping[m_endpoint_indices[b].component[cAlpha0]];
        if (!m_endpoint_indices[b].reference) {
          int sym = index - index0;
          hist[sym < 0 ? sym + codebook_size : sym]++;
        }
        index0 = index;
      }
      if (hasAlpha1) {
        uint index = remapping[m_endpoint_indices[b].component[cAlpha1]];
        if (!m_endpoint_indices[b].reference) {
          int sym = index - index1;
          hist[sym < 0 ? sym + codebook_size : sym]++;
        }
        index1 = index;
      }
    }
  }

  static_huffman_data_model dm;
  dm.init(true, codebook_size, hist.get_ptr(), 16);
  const uint8* code_sizes = dm.get_code_sizes();
  for (uint s = 0; s < codebook_size; s++)
    total_bits += hist[s] * code_sizes[s];

  symbol_codec codec;
  codec.start_encoding(64 * 1024);
  codec.encode_enable_simulation(true);
  codec.encode_transmit_static_huffman_data_model(dm, false);
  codec.stop_encoding(false);
  total_bits += codec.encode_get_total_bits_written();

  pParams->m_trial->total_bits = total_bits;

  crnlib_delete(pParams);
}

bool crn_comp::optimize_alpha_endpoint_codebook(crnlib::vector<uint>& remapping) {
  if (m_pParams->m_flags & cCRNCompFlagQuick) {
    remapping.resize(m_hvq.get_alpha_endpoint_vec().size());
    for (uint i = 0; i < m_hvq.get_alpha_endpoint_vec().size(); i++)
      remapping[i] = i;

    if (!pack_alpha_endpoints(m_packed_alpha_endpoints, remapping, 0))
      return false;

    return true;
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("----- Begin optimization of alpha endpoint codebook");
#endif

  const uint cMaxEndpointRemapIters = 3;
  optimize_alpha_endpoint_codebook_params::trial remapping_trial[cMaxEndpointRemapIters + 1];

  uint n = m_hvq.get_alpha_endpoint_codebook_size();
  hist_type xhist(n * n);
  bool hasAlpha0 = m_has_comp[cAlpha0], hasAlpha1 = m_has_comp[cAlpha1];
  for (uint b = 1; b < m_endpoint_indices.size(); b++) {
    if (!m_endpoint_indices[b].reference) {
      if (hasAlpha0) {
        update_hist(xhist, m_endpoint_indices[b - 1].alpha0, m_endpoint_indices[b].alpha0, n);
        update_hist(xhist, m_endpoint_indices[b].alpha0, m_endpoint_indices[b - 1].alpha0, n);
      }
      if (hasAlpha1) {
        update_hist(xhist, m_endpoint_indices[b - 1].alpha1, m_endpoint_indices[b].alpha1, n);
        update_hist(xhist, m_endpoint_indices[b].alpha1, m_endpoint_indices[b - 1].alpha1, n);
      }
    }
  }

  for (uint i = 0; i <= cMaxEndpointRemapIters; i++) {
    optimize_alpha_endpoint_codebook_params* pParams = crnlib_new<optimize_alpha_endpoint_codebook_params>();
    pParams->m_iter_index = i;
    pParams->m_max_iter_index = cMaxEndpointRemapIters;
    pParams->m_trial = remapping_trial + i;
    pParams->m_xhist = &xhist;
    m_task_pool.queue_object_task(this, &crn_comp::optimize_alpha_endpoint_codebook_task, 0, pParams);
  }
  m_task_pool.join();

  for (uint best_bits = UINT_MAX, i = 0; i <= cMaxEndpointRemapIters; i++) {
    if (remapping_trial[i].total_bits < best_bits) {
      m_packed_alpha_endpoints.swap(remapping_trial[i].packed_endpoints);
      remapping.swap(remapping_trial[i].remapping);
      best_bits = remapping_trial[i].total_bits;
    }
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    console::debug("End optimization of alpha endpoint codebook");
#endif

  return true;
}

bool crn_comp::optimize_alpha_selector_codebook(crnlib::vector<uint>& remapping) {
  if (m_pParams->m_flags & cCRNCompFlagQuick) {
    remapping.resize(m_hvq.get_alpha_selectors_vec().size());
    for (uint i = 0; i < m_hvq.get_alpha_selectors_vec().size(); i++)
      remapping[i] = i;
  } else {
    sort_selector_codebook(remapping, m_hvq.get_alpha_selectors_vec(), g_dxt5_to_linear);
  }
  return pack_selectors(m_packed_alpha_selectors, m_hvq.get_alpha_selectors_vec(), remapping, 7, g_dxt5_to_linear, 0);
}

bool crn_comp::pack_data_models() {
  symbol_codec codec;
  codec.start_encoding(1024 * 1024);

  if (!codec.encode_transmit_static_huffman_data_model(m_reference_encoding_dm, false))
    return false;

  for (uint i = 0; i < 2; i++) {
    if (m_endpoint_index_dm[i].get_total_syms()) {
      if (!codec.encode_transmit_static_huffman_data_model(m_endpoint_index_dm[i], false))
        return false;
    }

    if (m_selector_index_dm[i].get_total_syms()) {
      if (!codec.encode_transmit_static_huffman_data_model(m_selector_index_dm[i], false))
        return false;
    }
  }

  codec.stop_encoding(false);

  m_packed_data_models.swap(codec.get_encoding_buf());

  return true;
}

bool crn_comp::create_comp_data() {
  utils::zero_object(m_crn_header);

  m_crn_header.m_width = static_cast<uint16>(m_pParams->m_width);
  m_crn_header.m_height = static_cast<uint16>(m_pParams->m_height);
  m_crn_header.m_levels = static_cast<uint8>(m_pParams->m_levels);
  m_crn_header.m_faces = static_cast<uint8>(m_pParams->m_faces);
  m_crn_header.m_format = static_cast<uint8>(m_pParams->m_format);
  m_crn_header.m_userdata0 = m_pParams->m_userdata0;
  m_crn_header.m_userdata1 = m_pParams->m_userdata1;

  m_comp_data.clear();
  m_comp_data.reserve(2 * 1024 * 1024);
  append_vec(m_comp_data, &m_crn_header, sizeof(m_crn_header));
  // tack on the rest of the variable size m_level_ofs array
  m_comp_data.resize(m_comp_data.size() + sizeof(m_crn_header.m_level_ofs[0]) * (m_pParams->m_levels - 1));

  if (m_packed_color_endpoints.size()) {
    m_crn_header.m_color_endpoints.m_num = static_cast<uint16>(m_hvq.get_color_endpoint_codebook_size());
    m_crn_header.m_color_endpoints.m_size = m_packed_color_endpoints.size();
    m_crn_header.m_color_endpoints.m_ofs = m_comp_data.size();
    append_vec(m_comp_data, m_packed_color_endpoints);
  }

  if (m_packed_color_selectors.size()) {
    m_crn_header.m_color_selectors.m_num = static_cast<uint16>(m_hvq.get_color_selector_codebook_size());
    m_crn_header.m_color_selectors.m_size = m_packed_color_selectors.size();
    m_crn_header.m_color_selectors.m_ofs = m_comp_data.size();
    append_vec(m_comp_data, m_packed_color_selectors);
  }

  if (m_packed_alpha_endpoints.size()) {
    m_crn_header.m_alpha_endpoints.m_num = static_cast<uint16>(m_hvq.get_alpha_endpoint_codebook_size());
    m_crn_header.m_alpha_endpoints.m_size = m_packed_alpha_endpoints.size();
    m_crn_header.m_alpha_endpoints.m_ofs = m_comp_data.size();
    append_vec(m_comp_data, m_packed_alpha_endpoints);
  }

  if (m_packed_alpha_selectors.size()) {
    m_crn_header.m_alpha_selectors.m_num = static_cast<uint16>(m_hvq.get_alpha_selector_codebook_size());
    m_crn_header.m_alpha_selectors.m_size = m_packed_alpha_selectors.size();
    m_crn_header.m_alpha_selectors.m_ofs = m_comp_data.size();
    append_vec(m_comp_data, m_packed_alpha_selectors);
  }

  m_crn_header.m_tables_ofs = m_comp_data.size();
  m_crn_header.m_tables_size = m_packed_data_models.size();
  append_vec(m_comp_data, m_packed_data_models);

  uint level_ofs[cCRNMaxLevels];
  for (uint i = 0; i < m_mip_groups.size(); i++) {
    level_ofs[i] = m_comp_data.size();
    append_vec(m_comp_data, m_packed_chunks[i]);
  }

  crnd::crn_header& dst_header = *(crnd::crn_header*)&m_comp_data[0];
  // don't change the m_comp_data vector - or dst_header will be invalidated!

  memcpy(&dst_header, &m_crn_header, sizeof(dst_header));

  for (uint i = 0; i < m_mip_groups.size(); i++)
    dst_header.m_level_ofs[i] = level_ofs[i];

  const uint actual_header_size = sizeof(crnd::crn_header) + sizeof(dst_header.m_level_ofs[0]) * (m_mip_groups.size() - 1);

  dst_header.m_sig = crnd::crn_header::cCRNSigValue;

  dst_header.m_data_size = m_comp_data.size();
  dst_header.m_data_crc16 = crc16(&m_comp_data[actual_header_size], m_comp_data.size() - actual_header_size);

  dst_header.m_header_size = actual_header_size;
  dst_header.m_header_crc16 = crc16(&dst_header.m_data_size, actual_header_size - (uint)((uint8*)&dst_header.m_data_size - (uint8*)&dst_header));

  return true;
}

bool crn_comp::update_progress(uint phase_index, uint subphase_index, uint subphase_total) {
  if (!m_pParams->m_pProgress_func)
    return true;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_pParams->m_flags & cCRNCompFlagDebugging)
    return true;
#endif

  return (*m_pParams->m_pProgress_func)(phase_index, cTotalCompressionPhases, subphase_index, subphase_total, m_pParams->m_pProgress_func_data) != 0;
}

bool crn_comp::compress_internal() {
  if (!alias_images())
    return false;
  if (!quantize_images())
    return false;

  crnlib::vector<uint> endpoint_remap[2];
  crnlib::vector<uint> selector_remap[2];

  if (m_has_comp[cColor]) {
    if (!optimize_color_endpoint_codebook(endpoint_remap[0]))
      return false;
    if (!optimize_color_selector_codebook(selector_remap[0]))
      return false;
  }

  if (m_has_comp[cAlpha0]) {
    if (!optimize_alpha_endpoint_codebook(endpoint_remap[1]))
      return false;
    if (!optimize_alpha_selector_codebook(selector_remap[1]))
      return false;
  }

  m_chunk_encoding_hist.clear();
  for (uint i = 0; i < 2; i++) {
    m_endpoint_index_hist[i].clear();
    m_endpoint_index_dm[i].clear();
    m_selector_index_hist[i].clear();
    m_selector_index_dm[i].clear();
  }

  for (uint pass = 0; pass < 2; pass++) {
    for (uint mip_group = 0; mip_group < m_mip_groups.size(); mip_group++) {
      symbol_codec codec;
      codec.start_encoding(2 * 1024 * 1024);

      if (!pack_chunks(
          mip_group,
          !pass && !mip_group, pass ? &codec : NULL,
          m_has_comp[cColor] ? &endpoint_remap[0] : NULL, m_has_comp[cColor] ? &selector_remap[0] : NULL,
          m_has_comp[cAlpha0] ? &endpoint_remap[1] : NULL, m_has_comp[cAlpha0] ? &selector_remap[1] : NULL)) {
        return false;
      }

      codec.stop_encoding(false);

      if (pass)
        m_packed_chunks[mip_group].swap(codec.get_encoding_buf());
    }

    if (!pass) {
      m_reference_encoding_dm.init(true, m_chunk_encoding_hist, 16);

      for (uint i = 0; i < 2; i++) {
        if (m_endpoint_index_hist[i].size())
          m_endpoint_index_dm[i].init(true, m_endpoint_index_hist[i], 16);

        if (m_selector_index_hist[i].size())
          m_selector_index_dm[i].init(true, m_selector_index_hist[i], 16);
      }
    }
  }

  if (!pack_data_models())
    return false;

  if (!create_comp_data())
    return false;

  if (!update_progress(24, 1, 1))
    return false;

  if (m_pParams->m_flags & cCRNCompFlagDebugging) {
    crnlib_print_mem_stats();
  }

  return true;
}

bool crn_comp::compress_init(const crn_comp_params& params) {
  params;
  return true;
}

bool crn_comp::compress_pass(const crn_comp_params& params, float* pEffective_bitrate) {
  clear();

  if (pEffective_bitrate)
    *pEffective_bitrate = 0.0f;

  m_pParams = &params;

  if ((math::minimum(m_pParams->m_width, m_pParams->m_height) < 1) || (math::maximum(m_pParams->m_width, m_pParams->m_height) > cCRNMaxLevelResolution))
    return false;

  if (!m_task_pool.init(params.m_num_helper_threads))
    return false;

  bool status = compress_internal();

  m_task_pool.deinit();

  if ((status) && (pEffective_bitrate)) {
    uint total_pixels = 0;

    for (uint f = 0; f < m_pParams->m_faces; f++)
      for (uint l = 0; l < m_pParams->m_levels; l++)
        total_pixels += m_images[f][l].get_total_pixels();

    *pEffective_bitrate = (m_comp_data.size() * 8.0f) / total_pixels;
  }

  return status;
}

void crn_comp::compress_deinit() {
}

}  // namespace crnlib
