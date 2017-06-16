// File: crn_dxt_hc.cpp
// See Copyright Notice and license at the end of inc/crnlib.h
#include "crn_core.h"
#include "crn_dxt_hc.h"
#include "crn_image_utils.h"
#include "crn_console.h"
#include "crn_dxt_fast.h"

namespace crnlib {

typedef vec<6, float> vec6F;
typedef vec<16, float> vec16F;

static uint8 g_tile_map[8][2][2] = {
  {{ 0, 0 }, { 0, 0 }},
  {{ 0, 0 }, { 1, 1 }},
  {{ 0, 1 }, { 0, 1 }},
  {{ 0, 0 }, { 1, 2 }},
  {{ 1, 2 }, { 0, 0 }},
  {{ 0, 1 }, { 0, 2 }},
  {{ 1, 0 }, { 2, 0 }},
  {{ 0, 1 }, { 2, 3 }},
};

dxt_hc::dxt_hc()
  : m_num_blocks(0),
    m_has_color_blocks(false),
    m_num_alpha_blocks(0),
    m_main_thread_id(crn_get_current_thread_id()),
    m_canceled(false),
    m_pTask_pool(NULL),
    m_prev_phase_index(-1),
    m_prev_percentage_complete(-1) {
}

dxt_hc::~dxt_hc() {
}

void dxt_hc::clear() {
  m_blocks = 0;
  m_num_blocks = 0;
  m_num_alpha_blocks = 0;
  m_has_color_blocks = false;

  m_color_clusters.clear();
  m_alpha_clusters.clear();

  m_canceled = false;

  m_prev_phase_index = -1;
  m_prev_percentage_complete = -1;

  m_block_weights.clear();
  m_block_encodings.clear();
  for (uint c = 0; c < 3; c++)
    m_block_selectors[c].clear();
  m_color_selectors.clear();
  m_alpha_selectors.clear();
  m_color_selectors_used.clear();
  m_alpha_selectors_used.clear();
  m_tile_indices.clear();
  m_endpoint_indices.clear();
  m_selector_indices.clear();
  m_tiles.clear();
  m_num_tiles = 0;
}

bool dxt_hc::compress(
    color_quad_u8 (*blocks)[16],
    crnlib::vector<endpoint_indices_details>& endpoint_indices,
    crnlib::vector<selector_indices_details>& selector_indices,
    crnlib::vector<uint32>& color_endpoints,
    crnlib::vector<uint32>& alpha_endpoints,
    crnlib::vector<uint32>& color_selectors,
    crnlib::vector<uint64>& alpha_selectors,
    const params& p
  ) {
  clear();
  m_has_color_blocks = p.m_format == cDXT1 || p.m_format == cDXT5;
  m_num_alpha_blocks = p.m_format == cDXT5 || p.m_format == cDXT5A ? 1 : p.m_format == cDXN_XY || p.m_format == cDXN_YX ? 2 : 0;
  if (!m_has_color_blocks && !m_num_alpha_blocks)
    return false;
  m_blocks = blocks;
  m_main_thread_id = crn_get_current_thread_id();
  m_pTask_pool = p.m_pTask_pool;
  m_params = p;

  uint tile_derating[8] = {0, 1, 1, 2, 2, 2, 2, 3};
  for (uint level = 0; level < p.m_num_levels; level++) {
    float adaptive_tile_color_psnr_derating = p.m_adaptive_tile_color_psnr_derating;
    if (level && adaptive_tile_color_psnr_derating > .25f)
      adaptive_tile_color_psnr_derating = math::maximum(.25f, adaptive_tile_color_psnr_derating / powf(3.0f, static_cast<float>(level)));
    for (uint e = 0; e < 8; e++)
      m_color_derating[level][e] = math::lerp(0.0f, adaptive_tile_color_psnr_derating, tile_derating[e] / 3.0f);
  }
  for (uint e = 0; e < 8; e++)
    m_alpha_derating[e] = math::lerp(0.0f, m_params.m_adaptive_tile_alpha_psnr_derating, tile_derating[e] / 3.0f);
  for (uint i = 0; i < 256; i++)
    m_uint8_to_float[i] = i * 1.0f / 255.0f;

  m_num_blocks = m_params.m_num_blocks;
  m_block_weights.resize(m_num_blocks);
  m_block_encodings.resize(m_num_blocks);
  for (uint c = 0; c < 3; c++)
    m_block_selectors[c].resize(m_num_blocks);
  m_tile_indices.resize(m_num_blocks);
  m_endpoint_indices.resize(m_num_blocks);
  m_selector_indices.resize(m_num_blocks);
  m_tiles.resize(m_num_blocks);

  for (uint level = 0; level < p.m_num_levels; level++) {
    float weight = p.m_levels[level].m_weight;
    for (uint b = p.m_levels[level].m_first_block, bEnd = b + p.m_levels[level].m_num_blocks; b < bEnd; b++)
      m_block_weights[b] = weight;
  }

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::determine_tiles_task, i);
  m_pTask_pool->join();

  m_num_tiles = 0;
  for (uint t = 0; t < m_tiles.size(); t++) {
    if (m_tiles[t].pixels.size())
      m_num_tiles++;
  }

  if (m_has_color_blocks)
    determine_color_endpoints();

  if (m_num_alpha_blocks)
    determine_alpha_endpoints();

  if (m_has_color_blocks)
    create_color_selector_codebook();

  if (m_num_alpha_blocks)
    create_alpha_selector_codebook();

  color_endpoints.reserve(color_endpoints.size() + m_color_clusters.size());
  crnlib::vector<uint16> color_endpoints_remap(m_color_clusters.size());
  hash_map<uint32, uint> color_endpoints_map;
  for (uint i = 0; i < m_color_clusters.size(); i++) {
    if (m_color_clusters[i].pixels.size()) {
      uint32 endpoint = dxt1_block::pack_endpoints(m_color_clusters[i].first_endpoint, m_color_clusters[i].second_endpoint);
      hash_map<uint32, uint>::insert_result insert_result = color_endpoints_map.insert(endpoint, color_endpoints.size());
      if (insert_result.second) {
        color_endpoints_remap[i] = color_endpoints.size();
        color_endpoints.push_back(endpoint);
      } else {
        color_endpoints_remap[i] = insert_result.first->second;
      }
    }
  }

  alpha_endpoints.reserve(alpha_endpoints.size() + m_alpha_clusters.size());
  crnlib::vector<uint16> alpha_endpoints_remap(m_alpha_clusters.size());
  hash_map<uint32, uint> alpha_endpoints_map;
  for (uint i = 0; i < m_alpha_clusters.size(); i++) {
    if (m_alpha_clusters[i].pixels.size()) {
      uint32 endpoint = dxt5_block::pack_endpoints(m_alpha_clusters[i].first_endpoint, m_alpha_clusters[i].second_endpoint);
      hash_map<uint32, uint>::insert_result insert_result = alpha_endpoints_map.insert(endpoint, alpha_endpoints.size());
      if (insert_result.second) {
        alpha_endpoints_remap[i] = alpha_endpoints.size();
        alpha_endpoints.push_back(endpoint);
      } else {
        alpha_endpoints_remap[i] = insert_result.first->second;
      }
    }
  }

  color_selectors.reserve(color_selectors.size() + m_color_selectors.size());
  crnlib::vector<uint16> color_selectors_remap(m_color_selectors.size());
  hash_map<uint32, uint> color_selectors_map;
  for (uint i = 0; i < m_color_selectors.size(); i++) {
    if (m_color_selectors_used[i]) {
      hash_map<uint32, uint>::insert_result insert_result = color_selectors_map.insert(m_color_selectors[i], color_selectors.size());
      if (insert_result.second) {
        color_selectors_remap[i] = color_selectors.size();
        color_selectors.push_back(m_color_selectors[i]);
      } else {
        color_selectors_remap[i] = insert_result.first->second;
      }
    }
  }

  alpha_selectors.reserve(alpha_selectors.size() + m_alpha_selectors.size());
  crnlib::vector<uint16> alpha_selectors_remap(m_alpha_selectors.size());
  hash_map<uint64, uint> alpha_selectors_map;
  for (uint i = 0; i < m_alpha_selectors.size(); i++) {
    if (m_alpha_selectors_used[i]) {
      hash_map<uint64, uint>::insert_result insert_result = alpha_selectors_map.insert(m_alpha_selectors[i], alpha_selectors.size());
      if (insert_result.second) {
        alpha_selectors_remap[i] = alpha_selectors.size();
        alpha_selectors.push_back(m_alpha_selectors[i]);
      } else {
        alpha_selectors_remap[i] = insert_result.first->second;
      }
    }
  }

  endpoint_indices.resize(m_num_blocks);
  selector_indices.resize(m_num_blocks);
  for (uint level = 0; level < p.m_num_levels; level++) {
    uint first_block = p.m_levels[level].m_first_block;
    uint end_block = first_block + p.m_levels[level].m_num_blocks;
    uint block_width = p.m_levels[level].m_block_width;
    for (uint by = 0, b = first_block; b < end_block; by++) {
      for (uint bx = 0; bx < block_width; bx++, b++) {
        bool top_match = by != 0;
        bool left_match = top_match || bx;
        for (uint c = m_has_color_blocks ? 0 : cAlpha0; c < cAlpha0 + m_num_alpha_blocks; c++) {
          uint16 endpoint_index = (c ? alpha_endpoints_remap : color_endpoints_remap)[m_endpoint_indices[b].component[c]];
          left_match = left_match && endpoint_index == endpoint_indices[b - 1].component[c];
          top_match = top_match && endpoint_index == endpoint_indices[b - block_width].component[c];
          endpoint_indices[b].component[c] = endpoint_index;
          uint16 selector_index = (c ? alpha_selectors_remap : color_selectors_remap)[m_selector_indices[b].component[c]];
          selector_indices[b].component[c] = selector_index;
        }
        endpoint_indices[b].reference = left_match ? 1 : top_match ? 2 : 0;
      }
    }
  }

  m_pTask_pool = NULL;
  return true;
}

void dxt_hc::determine_tiles_task(uint64 data, void* pData_ptr) {
  uint num_tasks = m_pTask_pool->get_num_threads() + 1;
  uint offsets[9] = {0, 16, 32, 48, 0, 32, 64, 96, 64};
  uint8 tiles[8][4] = {{8}, {6, 7}, {4, 5}, {6, 1, 3}, {7, 0, 2}, {4, 2, 3}, {5, 0, 1}, {0, 2, 1, 3}};
  color_quad_u8 tilePixels[128];
  uint8 selectors[64];
  uint tile_error[3][9];
  uint total_error[3][8];
  tree_clusterizer<vec3F> color_palettizer;
  tree_clusterizer<vec1F> alpha_palettizer;

  for (uint level = 0; level < m_params.m_num_levels; level++) {
    float weight = m_params.m_levels[level].m_weight;
    uint width = m_params.m_levels[level].m_block_width;
    uint height = m_params.m_levels[level].m_num_blocks / width;
    uint faceHeight = height / m_params.m_num_faces;
    uint h = height * data / num_tasks & ~1;
    uint hEnd = height * (data + 1) / num_tasks & ~1;
    uint hFace = h % faceHeight;
    uint b = m_params.m_levels[level].m_first_block + h * width;

    for (; h < hEnd; h += 2, hFace += 2, b += width) {
      uint tile_offset = b;
      uint tile_offset_delta = 4;
      if (hFace == faceHeight) {
        hFace = 0;
      } else if (hFace & 2) {
        tile_offset_delta = -4;
        tile_offset += (width << 1) + tile_offset_delta;
      }
      for (uint bNext = b + width; b < bNext; b += 2, tile_offset += tile_offset_delta) {
        for (int t = 0; t < 64; t += 16)
          memcpy(tilePixels + t, m_blocks[b + (t & 16 ? width : 0) + (t & 32 ? 1 : 0)], 64);
        for (int t = 0; t < 64; t += 4)
          memcpy(tilePixels + 64 + t, m_blocks[b + (t & 32 ? width : 0) + (t & 4 ? 1 : 0)] + (t >> 1 & 12), 16);

        for (uint t = 0; t < 9; t++) {
          color_quad_u8* pixels = tilePixels + offsets[t];
          uint size = 16 << (t >> 2);
          if (m_has_color_blocks) {
            uint low16, high16;
            dxt_fast::compress_color_block(size, pixels, low16, high16, selectors);
            color_quad_u8 block_colors[4];
            dxt1_block::get_block_colors4(block_colors, low16, high16);
            uint error = 0;
            for (uint p = 0; p < size; p++) {
              for (uint8 c = 0; c < 3; c++) {
                uint delta = pixels[p][c] - block_colors[selectors[p]][c];
                error += delta * delta;
              }
            }
            tile_error[cColor][t] = error;
          }
          for (uint a = 0; a < m_num_alpha_blocks; a++) {
            uint8 component = m_params.m_alpha_component_indices[a];
            dxt5_endpoint_optimizer optimizer;
            dxt5_endpoint_optimizer::params params;
            dxt5_endpoint_optimizer::results results;
            params.m_pPixels = pixels;
            params.m_num_pixels = size;
            params.m_comp_index = component;
            params.m_use_both_block_types = false;
            params.m_quality = cCRNDXTQualityNormal;
            results.m_pSelectors = selectors;
            optimizer.compute(params, results);
            uint block_values[cDXT5SelectorValues];
            dxt5_block::get_block_values8(block_values, results.m_first_endpoint, results.m_second_endpoint);
            tile_error[cAlpha0 + a][t] = results.m_error;
          }
        }

        for (uint8 c = m_has_color_blocks ? 0 : cAlpha0; c < cAlpha0 + m_num_alpha_blocks; c++) {
          for (uint8 e = 0; e < 8; e++) {
            total_error[c][e] = 0;
            for (uint8 t = 0, s = e + 1; s; s >>= 1, t++)
              total_error[c][e] += tile_error[c][tiles[e][t]];
          }
        }

        float best_quality = 0.0f;
        uint best_encoding = 0;
        for (uint e = 0; e < 8; e++) {
          float quality = 0;
          if (m_has_color_blocks) {
            double peakSNR = total_error[cColor][e] ? log10(255.0f / sqrt(total_error[cColor][e] / 192.0)) * 20.0f : 999999.0f;
            quality = (float)math::maximum<double>(peakSNR - m_color_derating[level][e], 0.0f);
            if (m_num_alpha_blocks)
              quality *= m_params.m_adaptive_tile_color_alpha_weighting_ratio;
          }
          for (uint a = 0; a < m_num_alpha_blocks; a++) {
            double peakSNR = total_error[cAlpha0 + a][e] ? log10(255.0f / sqrt(total_error[cAlpha0 + a][e] / 64.0)) * 20.0f : 999999.0f;
            quality += (float)math::maximum<double>(peakSNR - m_alpha_derating[e], 0.0f);
          }
          if (quality > best_quality) {
            best_quality = quality;
            best_encoding = e;
          }
        }
    
        for (uint tile_index = 0, s = best_encoding + 1; s; s >>= 1, tile_index++) {
          tile_details& tile = m_tiles[tile_offset | tile_index];
          uint t = tiles[best_encoding][tile_index];
          tile.pixels.append(tilePixels + offsets[t], 16 << (t >> 2));
          tile.weight = weight;

          if (m_has_color_blocks) {
            color_palettizer.clear();
            for (uint p = 0; p < tile.pixels.size(); p++) {
              const color_quad_u8& pixel = tile.pixels[p];
              vec3F v(m_uint8_to_float[pixel[0]], m_uint8_to_float[pixel[1]], m_uint8_to_float[pixel[2]]);
              color_palettizer.add_training_vec(m_params.m_perceptual ? vec3F(v[0] * 0.5f, v[1], v[2] * 0.25f): v, 1);
            }
            color_palettizer.generate_codebook(2);
            bool single = color_palettizer.get_codebook_size() == 1;
            bool reorder = !single && color_palettizer.get_codebook_entry(0).length() > color_palettizer.get_codebook_entry(1).length();
            for (uint t = 0, i = 0; i < 2; i++) {
              vec3F v = color_palettizer.get_codebook_entry(single ? 0 : reorder ? 1 - i : i);
              for (uint c = 0; c < 3; c++, t++)
                tile.color_endpoint[t] = v[c];
            }
          }

          for (uint a = 0; a < m_num_alpha_blocks; a++) {
            alpha_palettizer.clear();
            for (uint c = m_params.m_alpha_component_indices[a], p = 0; p < tile.pixels.size(); p++)
              alpha_palettizer.add_training_vec(vec1F(m_uint8_to_float[tile.pixels[p][c]]), 1);
            alpha_palettizer.generate_codebook(2);
            float v[2] = {alpha_palettizer.get_codebook_entry(0)[0], alpha_palettizer.get_codebook_entry(alpha_palettizer.get_codebook_size() - 1)[0]};
            tile.alpha_endpoints[a][0] = math::minimum(v[0], v[1]);
            tile.alpha_endpoints[a][1] = math::maximum(v[0], v[1]);
          }
        }

        for (uint by = 0; by < 2; by++) {
          for (uint bx = 0; bx < 2; bx++) {
            m_block_encodings[b + (by ? width : 0) + bx] = best_encoding;
            m_tile_indices[b + (by ? width : 0) + bx] = tile_offset | g_tile_map[best_encoding][by][bx];
          }
        }

      }
    }
  }
}

void dxt_hc::determine_color_endpoint_codebook_task(uint64 data, void* pData_ptr) {
  pData_ptr;
  const uint thread_index = static_cast<uint>(data);

  if (!m_has_color_blocks)
    return;

  uint total_empty_clusters = 0;
  for (uint cluster_index = 0; cluster_index < m_color_clusters.size(); cluster_index++) {
    if (m_canceled)
      return;

    if ((crn_get_current_thread_id() == m_main_thread_id) && ((cluster_index & 63) == 0)) {
      if (!update_progress(3, cluster_index, m_color_clusters.size()))
        return;
    }

    if (m_pTask_pool->get_num_threads()) {
      if ((cluster_index % (m_pTask_pool->get_num_threads() + 1)) != thread_index)
        continue;
    }

    color_cluster& cluster = m_color_clusters[cluster_index];
    if (cluster.pixels.empty())
      continue;

    crnlib::vector<uint8> selectors(cluster.pixels.size());

    dxt1_endpoint_optimizer::params params;
    params.m_block_index = cluster_index;
    params.m_pPixels = cluster.pixels.get_ptr();
    params.m_num_pixels = cluster.pixels.size();
    params.m_pixels_have_alpha = false;
    params.m_use_alpha_blocks = false;
    params.m_perceptual = m_params.m_perceptual;
    params.m_quality = cCRNDXTQualityUber;
    params.m_endpoint_caching = false;

    dxt1_endpoint_optimizer::results results;
    results.m_pSelectors = selectors.get_ptr();

    dxt1_endpoint_optimizer optimizer;
    optimizer.compute(params, results);
    cluster.first_endpoint = results.m_low_color;
    cluster.second_endpoint = results.m_high_color;
    color_quad_u8 block_values[4], color_values[4];
    dxt1_block::get_block_colors4(block_values, cluster.first_endpoint, cluster.second_endpoint);
    for (uint i = 0; i < 4; i++)
      color_values[i] = cluster.color_values[i] = block_values[g_dxt1_from_linear[i]];
    for (uint c = 0; results.m_alternate_rounding && c < 3; c++) {
      color_values[1].c[c] = ((color_values[0].c[c] << 1) + color_values[3].c[c] + 1) / 3;
      color_values[2].c[c] = ((color_values[3].c[c] << 1) + color_values[0].c[c] + 1) / 3;
    }

    uint endpoint_weight = color::color_distance(m_params.m_perceptual, color_values[0], color_values[3], false) / 2000;
    float encoding_weight[8];
    for (uint i = 0; i < 8; i++)
      encoding_weight[i] = math::lerp(1.15f, 1.0f, i / 7.0f);

    crnlib::vector<uint>& blocks = cluster.blocks[cColor];
    for (uint i = 0; i < blocks.size(); i++) {
      uint b = blocks[i];
      uint weight = (uint)(math::clamp<uint>(endpoint_weight * m_block_weights[b], 1, 2048) * encoding_weight[m_block_encodings[b]]);
      uint32 selector = 0;
      for (uint sh = 0, p = 0; p < 16; p++, sh += 2) {
        uint error_best = cUINT32_MAX;
        uint8 s_best = 0;
        for (uint8 t = 0; t < 4; t++) {
          uint8 s = results.m_reordered ? 3 - g_dxt1_to_linear[t] : g_dxt1_to_linear[t];
          uint error = color::color_distance(m_params.m_perceptual, (color_quad_u8&)m_blocks[b][p], color_values[s], false);
          if (error < error_best) {
            s_best = s;
            error_best = error;
          }
        }
        selector |= s_best << sh;
      }
      m_block_selectors[cColor][b] = selector | (uint64)weight << 32;
    }

    dxt_endpoint_refiner refiner;
    dxt_endpoint_refiner::params refinerParams;
    dxt_endpoint_refiner::results refinerResults;
    refinerParams.m_perceptual = m_params.m_perceptual;
    refinerParams.m_pSelectors = selectors.get_ptr();
    refinerParams.m_pPixels = cluster.pixels.get_ptr();
    refinerParams.m_num_pixels = cluster.pixels.size();
    refinerParams.m_dxt1_selectors = true;
    refinerParams.m_error_to_beat = results.m_error;
    refinerParams.m_block_index = cluster_index;
    if (refiner.refine(refinerParams, refinerResults)) {
      cluster.first_endpoint = refinerResults.m_low_color;
      cluster.second_endpoint = refinerResults.m_high_color;
    }
  }
}

void dxt_hc::determine_color_endpoint_clusters_task(uint64 data, void* pData_ptr) {
  tree_clusterizer<vec6F>* vq = (tree_clusterizer<vec6F>*)pData_ptr;
  uint num_tasks = m_pTask_pool->get_num_threads() + 1;
  for (uint t = m_tiles.size() * data / num_tasks, tEnd = m_tiles.size() * (data + 1) / num_tasks; t < tEnd; t++) {
    if (m_tiles[t].pixels.size())
      m_tiles[t].cluster_indices[cColor] = vq->find_best_codebook_entry_fs(m_tiles[t].color_endpoint);
  }
}

void dxt_hc::determine_color_endpoints() {
  tree_clusterizer<vec6F> vq;
  for (uint t = 0; t < m_tiles.size(); t++) {
    if (m_tiles[t].pixels.size())
      vq.add_training_vec(m_tiles[t].color_endpoint, (uint)(m_tiles[t].pixels.size() * m_tiles[t].weight));
  }

  vq.generate_codebook(math::minimum<uint>(m_num_tiles, m_params.m_color_endpoint_codebook_size));
  m_color_clusters.resize(vq.get_codebook_size());

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::determine_color_endpoint_clusters_task, i, &vq);
  m_pTask_pool->join();

  for (uint t = 0; t < m_num_blocks; t++) {
    if (m_tiles[t].pixels.size())
      m_color_clusters[m_tiles[t].cluster_indices[cColor]].pixels.append(m_tiles[t].pixels);
  }

  for (uint b = 0; b < m_num_blocks; b++) {
    uint cluster_index = m_tiles[m_tile_indices[b]].cluster_indices[cColor];
    m_endpoint_indices[b].component[cColor] = cluster_index;
    m_color_clusters[cluster_index].blocks[cColor].push_back(b);
  }

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::determine_color_endpoint_codebook_task, i, NULL);
  m_pTask_pool->join();
}

void dxt_hc::determine_alpha_endpoint_codebook_task(uint64 data, void* pData_ptr) {
  pData_ptr;
  const uint thread_index = static_cast<uint>(data);

  for (uint cluster_index = 0; cluster_index < m_alpha_clusters.size(); cluster_index++) {
    if (m_canceled)
      return;

    if ((crn_get_current_thread_id() == m_main_thread_id) && ((cluster_index & 63) == 0)) {
      if (!update_progress(8, cluster_index, m_alpha_clusters.size()))
        return;
    }

    if (m_pTask_pool->get_num_threads()) {
      if ((cluster_index % (m_pTask_pool->get_num_threads() + 1)) != thread_index)
        continue;
    }

    alpha_cluster& cluster = m_alpha_clusters[cluster_index];
    if (cluster.pixels.empty())
      continue;

    crnlib::vector<uint8> selectors(cluster.pixels.size());

    dxt5_endpoint_optimizer::params params;
    params.m_pPixels = cluster.pixels.get_ptr();
    params.m_num_pixels = cluster.pixels.size();
    params.m_comp_index = 0;
    params.m_quality = cCRNDXTQualityUber;
    params.m_use_both_block_types = false;

    dxt5_endpoint_optimizer::results results;
    results.m_pSelectors = selectors.get_ptr();

    dxt5_endpoint_optimizer optimizer;
    optimizer.compute(params, results);
    cluster.first_endpoint = results.m_first_endpoint;
    cluster.second_endpoint = results.m_second_endpoint;
    uint block_values[8], alpha_values[8];
    dxt5_block::get_block_values(block_values, cluster.first_endpoint, cluster.second_endpoint);
    for (uint i = 0; i < 8; i++)
      alpha_values[i] = cluster.alpha_values[i] = block_values[g_dxt5_from_linear[i]];

    int delta = cluster.second_endpoint - cluster.first_endpoint;
    uint encoding_weight[8];
    for (uint endpoint_weight = math::clamp<uint>(delta * delta >> 3, 1, 2048), i = 0; i < 8; i++)
      encoding_weight[i] = (uint)(endpoint_weight * math::lerp(1.15f, 1.0f, i / 7.0f));

    for (uint a = 0; a < m_num_alpha_blocks; a++) {
      uint component_index = m_params.m_alpha_component_indices[a];
      crnlib::vector<uint>& blocks = cluster.blocks[cAlpha0 + a];
      for (uint i = 0; i < blocks.size(); i++) {
        uint b = blocks[i];
        uint weight = encoding_weight[m_block_encodings[b]];
        uint64 selector = 0;
        for (uint sh = 0, p = 0; p < 16; p++, sh += 3) {
          uint error_best = cUINT32_MAX;
          uint8 s_best = 0;
          for (uint8 t = 0; t < 8; t++) {
            uint8 s = results.m_reordered ? 7 - g_dxt5_to_linear[t] : g_dxt5_to_linear[t];
            int delta = m_blocks[b][p][component_index] - alpha_values[s];
            uint error = delta >= 0 ? delta : -delta;
            if (error < error_best) {
              s_best = s;
              error_best = error;
            }
          }
          selector |= (uint64)s_best << sh;
        }
        m_block_selectors[cAlpha0 + a][b] = selector | (uint64)weight << 48;
      }
    }

    dxt_endpoint_refiner refiner;
    dxt_endpoint_refiner::params refinerParams;
    dxt_endpoint_refiner::results refinerResults;
    refinerParams.m_perceptual = m_params.m_perceptual;
    refinerParams.m_pSelectors = selectors.get_ptr();
    refinerParams.m_pPixels = cluster.pixels.get_ptr();
    refinerParams.m_num_pixels = cluster.pixels.size();
    refinerParams.m_dxt1_selectors = false;
    refinerParams.m_error_to_beat = results.m_error;
    refinerParams.m_block_index = cluster_index;
    cluster.refined_alpha = refiner.refine(refinerParams, refinerResults);
    if (cluster.refined_alpha) {
      cluster.first_endpoint = refinerResults.m_low_color;
      cluster.second_endpoint = refinerResults.m_high_color;
      dxt5_block::get_block_values(block_values, cluster.first_endpoint, cluster.second_endpoint);
      for (uint i = 0; i < 8; i++)
        cluster.refined_alpha_values[i] = block_values[g_dxt5_from_linear[i]];
    } else {
      memcpy(cluster.refined_alpha_values, cluster.alpha_values, sizeof(cluster.refined_alpha_values));
    }
  }
}

void dxt_hc::determine_alpha_endpoint_clusters_task(uint64 data, void* pData_ptr) {
  tree_clusterizer<vec2F>* vq = (tree_clusterizer<vec2F>*)pData_ptr;
  uint num_tasks = m_pTask_pool->get_num_threads() + 1;
  for (uint t = m_tiles.size() * data / num_tasks, tEnd = m_tiles.size() * (data + 1) / num_tasks; t < tEnd; t++) {
    if (m_tiles[t].pixels.size()) {
      for (uint a = 0; a < m_num_alpha_blocks; a++)
        m_tiles[t].cluster_indices[cAlpha0 + a] = vq->find_best_codebook_entry_fs(m_tiles[t].alpha_endpoints[a]);
    }
  }
}

void dxt_hc::determine_alpha_endpoints() {
  tree_clusterizer<vec2F> vq;
  for (uint a = 0; a < m_num_alpha_blocks; a++) {
    for (uint t = 0; t < m_tiles.size(); t++) {
      if (m_tiles[t].pixels.size())
        vq.add_training_vec(m_tiles[t].alpha_endpoints[a], m_tiles[t].pixels.size());
    }
  }

  vq.generate_codebook(math::minimum<uint>(m_num_tiles, m_params.m_alpha_endpoint_codebook_size));
  m_alpha_clusters.resize(vq.get_codebook_size());

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::determine_alpha_endpoint_clusters_task, i, &vq);
  m_pTask_pool->join();

  for (uint a = 0; a < m_num_alpha_blocks; a++) {
    uint component_index = m_params.m_alpha_component_indices[a];
    for (uint t = 0; t < m_num_blocks; t++) {
      crnlib::vector<color_quad_u8>& source = m_tiles[t].pixels;
      if (source.size()) {
        crnlib::vector<color_quad_u8>& destination = m_alpha_clusters[m_tiles[t].cluster_indices[cAlpha0 + a]].pixels;
        for (uint p = 0; p < source.size(); p++)
          destination.push_back(color_quad_u8(source[p][component_index]));
      }
    }
  }

  for (uint b = 0; b < m_num_blocks; b++) {
    for (uint a = 0; a < m_num_alpha_blocks; a++) {
      uint cluster_index = m_tiles[m_tile_indices[b]].cluster_indices[cAlpha0 + a];
      m_endpoint_indices[b].component[cAlpha0 + a] = cluster_index;
      m_alpha_clusters[cluster_index].blocks[cAlpha0 + a].push_back(b);
    }
  }

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::determine_alpha_endpoint_codebook_task, i, NULL);
  m_pTask_pool->join();
}

struct color_selector_details {
  color_selector_details() { utils::zero_object(*this); }
  uint error[16][4];
  bool used;
};

void dxt_hc::create_color_selector_codebook_task(uint64 data, void* pData_ptr) {
  crnlib::vector<color_selector_details>& selector_details = *static_cast<crnlib::vector<color_selector_details>*>(pData_ptr);
  uint num_tasks = m_pTask_pool->get_num_threads() + 1;
  uint E2[16][4];
  uint E4[8][16];
  uint E8[4][256];
  for (uint b = m_num_blocks * data / num_tasks, bEnd = m_num_blocks * (data + 1) / num_tasks; b < bEnd; b++) {
    color_cluster& cluster = m_color_clusters[m_endpoint_indices[b].color];
    color_quad_u8* endpoint_colors = cluster.color_values;
    for (uint p = 0; p < 16; p++) {
      for (uint s = 0; s < 4; s++)
        E2[p][s] = color::color_distance(m_params.m_perceptual, m_blocks[b][p], endpoint_colors[s], false);
    }
    for (uint p = 0; p < 8; p++) {
      for (uint s = 0; s < 16; s++)
        E4[p][s] = E2[p << 1][s & 3] + E2[p << 1 | 1][s >> 2];
    }
    for (uint p = 0; p < 4; p++) {
      for (uint s = 0; s < 256; s++)
        E8[p][s] = E4[p << 1][s & 15] + E4[p << 1 | 1][s >> 4];
    }
    uint best_index = 0;
    for (uint best_error = cUINT32_MAX, s = 0; s < m_color_selectors.size(); s++) {
      uint32 selector = m_color_selectors[s];
      uint error = E8[0][selector & 255] + E8[1][selector >> 8 & 255] + E8[2][selector >> 16 & 255] + E8[3][selector >> 24 & 255];
      if (error < best_error) {
        best_error = error;
        best_index = s;
      }
    }
    uint (&total_errors)[16][4] = selector_details[best_index].error;
    for (uint p = 0; p < 16; p++) {
      for (uint s = 0; s < 4; s++)
        total_errors[p][s] += E2[p][s];
    }
    selector_details[best_index].used = true;
    m_selector_indices[b].color = best_index;
  }
}

void dxt_hc::create_color_selector_codebook() {
  tree_clusterizer<vec16F> selector_vq;
  vec16F v;
  for (uint b = 0; b < m_num_blocks; b++) {
    uint64 selector = m_block_selectors[cColor][b];
    for (uint8 p = 0; p < 16; p++, selector >>= 2)
      v[p] = ((selector & 3) + 0.5f) * 0.25f;
    selector_vq.add_training_vec(v, selector);
  }
  selector_vq.generate_codebook(m_params.m_color_selector_codebook_size);
  m_color_selectors.resize(selector_vq.get_codebook_size());
  m_color_selectors_used.resize(selector_vq.get_codebook_size());
  for (uint i = 0; i < selector_vq.get_codebook_size(); i++) {
    const vec16F& v = selector_vq.get_codebook_entry(i);
    m_color_selectors[i] = 0;
    for (uint sh = 0, j = 0; j < 16; j++, sh += 2)
      m_color_selectors[i] |= (uint)(v[j] * 4.0f) << sh;
  }

  uint num_tasks = m_pTask_pool->get_num_threads() + 1;
  crnlib::vector<crnlib::vector<color_selector_details>> selector_details(num_tasks);
  for (uint t = 0; t < num_tasks; t++) {
    selector_details[t].resize(m_color_selectors.size());
    m_pTask_pool->queue_object_task(this, &dxt_hc::create_color_selector_codebook_task, t, &selector_details[t]);
  }
  m_pTask_pool->join();

  for (uint t = 1; t < num_tasks; t++) {
    for (uint i = 0; i < m_color_selectors.size(); i++) {
      for (uint8 p = 0; p < 16; p++) {
        for (uint8 s = 0; s < 4; s++)
          selector_details[0][i].error[p][s] += selector_details[t][i].error[p][s];
      }
      selector_details[0][i].used = selector_details[0][i].used || selector_details[t][i].used;
    }
  }

  for (uint i = 0; i < m_color_selectors.size(); i++) {
    m_color_selectors_used[i] = selector_details[0][i].used;
    uint (&errors)[16][4] = selector_details[0][i].error;
    m_color_selectors[i] = 0;
    for (uint sh = 0, p = 0; p < 16; p++, sh += 2) {
      uint* e = errors[p];
      uint8 s03 = e[3] < e[0] ? 3 : 0;
      uint8 s12 = e[2] < e[1] ? 2 : 1;
      m_color_selectors[i] |= (e[s12] < e[s03] ? s12 : s03) << sh;
    }
  }
}

struct alpha_selector_details {
  alpha_selector_details() { utils::zero_object(*this); }
  uint error[16][8];
  bool used;
};

void dxt_hc::create_alpha_selector_codebook_task(uint64 data, void* pData_ptr) {
  crnlib::vector<alpha_selector_details>& selector_details = *static_cast<crnlib::vector<alpha_selector_details>*>(pData_ptr);
  uint num_tasks = m_pTask_pool->get_num_threads() + 1;
  uint E3[16][8];
  uint E6[8][64];
  for (uint b = m_num_blocks * data / num_tasks, bEnd = m_num_blocks * (data + 1) / num_tasks; b < bEnd; b++) {
    for (uint c = cAlpha0; c < cAlpha0 + m_num_alpha_blocks; c++) {
      const uint alpha_pixel_comp = m_params.m_alpha_component_indices[c - cAlpha0];
      alpha_cluster& cluster = m_alpha_clusters[m_endpoint_indices[b].component[c]];
      uint* block_values = cluster.alpha_values;
      for (uint p = 0; p < 16; p++) {
        for (uint s = 0; s < 8; s++) {
          int delta = m_blocks[b][p][alpha_pixel_comp] - block_values[s];
          E3[p][s] = delta * delta;
        }
      }
      for (uint p = 0; p < 8; p++) {
        for (uint s = 0; s < 64; s++)
          E6[p][s] = E3[p << 1][s & 7] + E3[p << 1 | 1][s >> 3];
      }
      uint best_index = 0;
      for (uint best_error = cUINT32_MAX, s = 0; s < m_alpha_selectors.size(); s++) {
        uint64 selector = m_alpha_selectors[s];
        uint error = E6[0][selector & 63];
        error += E6[1][selector >> 6 & 63];
        error += E6[2][selector >> 12 & 63];
        error += E6[3][selector >> 18 & 63];
        error += E6[4][selector >> 24 & 63];
        error += E6[5][selector >> 30 & 63];
        error += E6[6][selector >> 36 & 63];
        error += E6[7][selector >> 42 & 63];
        if (error < best_error) {
          best_error = error;
          best_index = s;
        }
      }
      if (cluster.refined_alpha) {
        block_values = cluster.refined_alpha_values;
        for (uint p = 0; p < 16; p++) {
          for (uint s = 0; s < 8; s++) {
            int delta = m_blocks[b][p][alpha_pixel_comp] - block_values[s];
            E3[p][s] = delta * delta;
          }
        }
      }
      uint (&total_errors)[16][8] = selector_details[best_index].error;
      for (uint p = 0; p < 16; p++) {
        for (uint s = 0; s < 8; s++)
         total_errors[p][s] += E3[p][s];
      }
      selector_details[best_index].used = true;
      m_selector_indices[b].component[c] = best_index;
    }
  }
}

void dxt_hc::create_alpha_selector_codebook() {
  tree_clusterizer<vec16F> selector_vq;
  vec16F v;
  for (uint c = cAlpha0; c < cAlpha0 + m_num_alpha_blocks; c++) {
    for (uint b = 0; b < m_num_blocks; b++) {
      uint64 selector = m_block_selectors[c][b];
      for (uint8 p = 0; p < 16; p++, selector >>= 3)
        v[p] = ((selector & 7) + 0.5f) * 0.125f;
      selector_vq.add_training_vec(v, selector);
    }
  }
  selector_vq.generate_codebook(m_params.m_alpha_selector_codebook_size);
  m_alpha_selectors.resize(selector_vq.get_codebook_size());
  m_alpha_selectors_used.resize(selector_vq.get_codebook_size());
  for (uint i = 0; i < selector_vq.get_codebook_size(); i++) {
    const vec16F& v = selector_vq.get_codebook_entry(i);
    m_alpha_selectors[i] = 0;
    for (uint sh = 0, j = 0; j < 16; j++, sh += 3)
      m_alpha_selectors[i] |= (uint64)(v[j] * 8.0f) << sh;
  }

  uint num_tasks = m_pTask_pool->get_num_threads() + 1;
  crnlib::vector<crnlib::vector<alpha_selector_details>> selector_details(num_tasks);
  for (uint t = 0; t < num_tasks; t++) {
    selector_details[t].resize(m_alpha_selectors.size());
    m_pTask_pool->queue_object_task(this, &dxt_hc::create_alpha_selector_codebook_task, t, &selector_details[t]);
  }
  m_pTask_pool->join();

  for (uint t = 1; t < num_tasks; t++) {
    for (uint i = 0; i < m_alpha_selectors.size(); i++) {
      for (uint8 p = 0; p < 16; p++) {
        for (uint8 s = 0; s < 8; s++)
          selector_details[0][i].error[p][s] += selector_details[t][i].error[p][s];
      }
      selector_details[0][i].used = selector_details[0][i].used || selector_details[t][i].used;
    }
  }

  for (uint i = 0; i < m_alpha_selectors.size(); i++) {
    m_alpha_selectors_used[i] = selector_details[0][i].used;
    uint (&errors)[16][8] = selector_details[0][i].error;
    m_alpha_selectors[i] = 0;
    for (uint sh = 0, p = 0; p < 16; p++, sh += 3) {
      uint* e = errors[p];
      uint8 s07 = e[7] < e[0] ? 7 : 0;
      uint8 s12 = e[2] < e[1] ? 2 : 1;
      uint8 s34 = e[4] < e[3] ? 4 : 3;
      uint8 s56 = e[6] < e[5] ? 6 : 5;
      uint8 s02 = e[s12] < e[s07] ? s12 : s07;
      uint8 s36 = e[s56] < e[s34] ? s56 : s34;
      m_alpha_selectors[i] |= (uint64)(e[s36] < e[s02] ? s36 : s02) << sh;
    }
  }
}

bool dxt_hc::update_progress(uint phase_index, uint subphase_index, uint subphase_total) {
  CRNLIB_ASSERT(crn_get_current_thread_id() == m_main_thread_id);

  if (!m_params.m_pProgress_func)
    return true;

  const int percentage_complete = (subphase_total > 1) ? ((100 * subphase_index) / (subphase_total - 1)) : 100;
  if (((int)phase_index == m_prev_phase_index) && (m_prev_percentage_complete == percentage_complete))
    return !m_canceled;

  m_prev_percentage_complete = percentage_complete;

  bool status = (*m_params.m_pProgress_func)(phase_index, cTotalCompressionPhases, subphase_index, subphase_total, m_params.m_pProgress_func_data) != 0;
  if (!status) {
    m_canceled = true;
    return false;
  }

  return true;
}

}  // namespace crnlib
