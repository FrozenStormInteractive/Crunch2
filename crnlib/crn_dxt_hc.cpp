// File: crn_dxt_hc.cpp
// See Copyright Notice and license at the end of inc/crnlib.h
#include "crn_core.h"
#include "crn_dxt_hc.h"
#include "crn_image_utils.h"
#include "crn_console.h"
#include "crn_dxt_fast.h"

#define CRNLIB_ENABLE_DEBUG_MESSAGES 0

namespace crnlib {

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
    m_num_alpha_blocks(0),
    m_has_color_blocks(false),
    m_has_alpha0_blocks(false),
    m_has_alpha1_blocks(false),
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
  m_has_alpha0_blocks = false;
  m_has_alpha1_blocks = false;

  m_color_clusters.clear();
  m_alpha_clusters.clear();
  m_alpha_selectors_vec.clear();
  m_color_selectors_vec.clear();

  m_color_endpoints.clear();
  m_alpha_endpoints.clear();

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

bool dxt_hc::compress(const params& p, task_pool& task_pool) {
  clear();
  m_params = p;
  m_main_thread_id = crn_get_current_thread_id();
  m_pTask_pool = &task_pool;

  switch (m_params.m_format) {
    case cDXT1: {
      m_has_color_blocks = true;
      break;
    }
    case cDXT5: {
      m_has_color_blocks = true;
      m_has_alpha0_blocks = true;
      m_num_alpha_blocks = 1;
      break;
    }
    case cDXT5A: {
      m_has_alpha0_blocks = true;
      m_num_alpha_blocks = 1;
      break;
    }
    case cDXN_XY:
    case cDXN_YX: {
      m_has_alpha0_blocks = true;
      m_has_alpha1_blocks = true;
      m_num_alpha_blocks = 2;
      break;
    }
    default: {
      return false;
    }
  }

  for (uint level = 0; level < p.m_num_levels; level++) {
    float adaptive_tile_color_psnr_derating = p.m_adaptive_tile_color_psnr_derating;
    if (level && adaptive_tile_color_psnr_derating > .25f)
      adaptive_tile_color_psnr_derating = math::maximum(.25f, adaptive_tile_color_psnr_derating / powf(3.0f, static_cast<float>(level)));
    for (uint e = 0; e < 8; e++)
      m_color_derating[level][e] = math::lerp(0.0f, adaptive_tile_color_psnr_derating, (g_chunk_encodings[e].m_num_tiles - 1) / 3.0f);
  }
  for (uint e = 0; e < 8; e++)
    m_alpha_derating[e] = math::lerp(0.0f, m_params.m_adaptive_tile_alpha_psnr_derating, (g_chunk_encodings[e].m_num_tiles - 1) / 3.0f);

  m_blocks = m_params.m_blocks;
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

  crnlib::vector<uint16> color_endpoint_remap(m_color_clusters.size());
  m_color_endpoints.reserve(m_color_clusters.size());
  hash_map<uint, uint> color_clusters_map;
  for (uint i = 0; i < m_color_clusters.size(); i++) {
    if (m_color_clusters[i].m_pixels.size()) {
      uint endpoint = dxt1_block::pack_endpoints(m_color_clusters[i].m_refined_first_endpoint, m_color_clusters[i].m_refined_second_endpoint);
      hash_map<uint, uint>::insert_result insert_result = color_clusters_map.insert(endpoint, m_color_endpoints.size());
      if (insert_result.second) {
        color_endpoint_remap[i] = m_color_endpoints.size();
        m_color_endpoints.push_back(endpoint);
      } else {
        color_endpoint_remap[i] = insert_result.first->second;
      }
    }
  }

  crnlib::vector<uint16> color_selector_remap(m_color_selectors.size());
  m_color_selectors_vec.reserve(m_color_selectors.size());
  hash_map<uint32, uint> color_selector_map;
  for (uint i = 0; i < m_color_selectors.size(); i++) {
    if (m_color_selectors_used[i]) {
      hash_map<uint32, uint>::insert_result insert_result = color_selector_map.insert(m_color_selectors[i], m_color_selectors_vec.size());
      if (insert_result.second) {
        color_selector_remap[i] = m_color_selectors_vec.size();
        selectors selector_vec;
        for (uint32 selector = m_color_selectors[i], s = 0; s < 16; s++, selector >>= 2)
          selector_vec.set_by_index(s, selector & 3);
        m_color_selectors_vec.push_back(selector_vec);
      } else {
        color_selector_remap[i] = insert_result.first->second;
      }
    }
  }

  crnlib::vector<uint16> alpha_endpoint_remap(m_alpha_clusters.size());
  m_alpha_endpoints.reserve(m_alpha_clusters.size());
  hash_map<uint, uint> alpha_endpoints_map;
  for (uint i = 0; i < m_alpha_clusters.size(); i++) {
    if (m_alpha_clusters[i].m_pixels.size()) {
      uint endpoint = dxt5_block::pack_endpoints(m_alpha_clusters[i].m_refined_first_endpoint, m_alpha_clusters[i].m_refined_second_endpoint);
      hash_map<uint, uint>::insert_result insert_result = alpha_endpoints_map.insert(endpoint, m_alpha_endpoints.size());
      if (insert_result.second) {
        alpha_endpoint_remap[i] = m_alpha_endpoints.size();
        m_alpha_endpoints.push_back(endpoint);
      } else {
        alpha_endpoint_remap[i] = insert_result.first->second;
      }
    }
  }

  crnlib::vector<uint16> alpha_selector_remap(m_alpha_selectors.size());
  m_alpha_selectors_vec.reserve(m_alpha_selectors.size());
  hash_map<uint64, uint> alpha_selectors_map;
  for (uint i = 0; i < m_alpha_selectors.size(); i++) {
    if (m_alpha_selectors_used[i]) {
      hash_map<uint64, uint>::insert_result insert_result = alpha_selectors_map.insert(m_alpha_selectors[i], m_alpha_selectors_vec.size());
      if (insert_result.second) {
        alpha_selector_remap[i] = m_alpha_selectors_vec.size();
        selectors selector_vec;
        for (uint64 selector = m_alpha_selectors[i], s = 0; s < 16; s++, selector >>= 3)
          selector_vec.set_by_index(s, selector & 7);
        m_alpha_selectors_vec.push_back(selector_vec);
      } else {
        alpha_selector_remap[i] = insert_result.first->second;
      }
    }
  }

  crnlib::vector<endpoint_indices_details>& endpoint_indices = *m_params.m_endpoint_indices;
  crnlib::vector<selector_indices_details>& selector_indices = *m_params.m_selector_indices;
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
        for (uint c = m_has_color_blocks ? 0 : cAlpha0Blocks; c < cAlpha0Blocks + m_num_alpha_blocks; c++) {
          uint16 endpoint_index = (c ? alpha_endpoint_remap : color_endpoint_remap)[m_endpoint_indices[b].component[c]];
          left_match = left_match && endpoint_index == endpoint_indices[b - 1].component[c];
          top_match = top_match && endpoint_index == endpoint_indices[b - block_width].component[c];
          endpoint_indices[b].component[c] = endpoint_index;
          uint16 selector_index = (c ? alpha_selector_remap : color_selector_remap)[m_selector_indices[b].component[c]];
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
  color_quad_u8 chunkPixels[128];
  uint8 selectors[64];
  uint tile_error[3][9];
  uint total_error[3][8];

  for (uint level = 0; level < m_params.m_num_levels; level++) {
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
          memcpy(chunkPixels + t, m_blocks[b + (t & 16 ? width : 0) + (t & 32 ? 1 : 0)], 64);
        for (int t = 0; t < 64; t += 4)
          memcpy(chunkPixels + 64 + t, m_blocks[b + (t & 32 ? width : 0) + (t & 4 ? 1 : 0)] + (t >> 1 & 12), 16);

        for (uint t = 0; t < 9; t++) {
          color_quad_u8* pixels = chunkPixels + offsets[t];
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
            tile_error[cColorBlocks][t] = error;
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
            tile_error[cAlpha0Blocks + a][t] = results.m_error;
          }
        }

        for (uint8 c = m_has_color_blocks ? 0 : cAlpha0Blocks; c < cAlpha0Blocks + m_num_alpha_blocks; c++) {
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
            double peakSNR = total_error[cColorBlocks][e] ? log10(255.0f / sqrt(total_error[cColorBlocks][e] / 192.0)) * 20.0f : 999999.0f;
            quality = (float)math::maximum<double>(peakSNR - m_color_derating[level][e], 0.0f);
            if (m_num_alpha_blocks)
              quality *= m_params.m_adaptive_tile_color_alpha_weighting_ratio;
          }
          for (uint a = 0; a < m_num_alpha_blocks; a++) {
            double peakSNR = total_error[cAlpha0Blocks + a][e] ? log10(255.0f / sqrt(total_error[cAlpha0Blocks + a][e] / 64.0)) * 20.0f : 999999.0f;
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
          tile.pixels.append(chunkPixels + offsets[t], 16 << (t >> 2));
          tile.weight = m_block_weights[b];

          if (m_has_color_blocks) {
            tree_clusterizer<vec3F> palettizer;
            for (uint p = 0; p < tile.pixels.size(); p++) {
              const color_quad_u8& c = tile.pixels[p];
              vec3F v(c[0] * 1.0f / 255.0f, c[1] * 1.0f / 255.0f, c[2] * 1.0f / 255.0f);
              if (m_params.m_perceptual) {
                v[0] *= 0.5f;
                v[2] *= 0.25f;
              }
              palettizer.add_training_vec(v, 1);
            }
            palettizer.generate_codebook(2);
            vec3F v[2];
            utils::zero_object(v);
            for (uint i = 0; i < palettizer.get_codebook_size(); i++)
              v[i] = palettizer.get_codebook_entry(i);
            if (palettizer.get_codebook_size() == 1)
              v[1] = v[0];
            if (v[0].length() > v[1].length())
              utils::swap(v[0], v[1]);
            vec6F vv;
            for (uint i = 0; i < 2; i++) {
              vv[i * 3 + 0] = v[i][0];
              vv[i * 3 + 1] = v[i][1];
              vv[i * 3 + 2] = v[i][2];
            }
            tile.color_endpoint = vv;
          }

          for (uint a = 0; a < m_num_alpha_blocks; a++) {
            uint component_index = m_params.m_alpha_component_indices[a];
            tree_clusterizer<vec1F> palettizer;
            for (uint p = 0; p < tile.pixels.size(); p++) {
              vec1F v(tile.pixels[p][component_index] * 1.0f / 255.0f);
              palettizer.add_training_vec(v, 1);
            }
            palettizer.generate_codebook(2);
            vec1F v[2];
            utils::zero_object(v);
            for (uint i = 0; i < palettizer.get_codebook_size(); i++)
              v[i] = palettizer.get_codebook_entry(i);
            if (palettizer.get_codebook_size() == 1)
              v[1] = v[0];
            if (v[0] > v[1])
              utils::swap(v[0], v[1]);
            vec2F vv(v[0][0], v[1][0]);
            tile.alpha_endpoints[a] = vv;
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

    endpoint_cluster& cluster = m_color_clusters[cluster_index];
    if (cluster.m_pixels.empty())
      continue;

    crnlib::vector<uint8> selectors(cluster.m_pixels.size());

    dxt1_endpoint_optimizer::params params;
    params.m_block_index = cluster_index;
    params.m_pPixels = cluster.m_pixels.get_ptr();
    params.m_num_pixels = cluster.m_pixels.size();
    params.m_pixels_have_alpha = false;
    params.m_use_alpha_blocks = false;
    params.m_perceptual = m_params.m_perceptual;
    params.m_quality = cCRNDXTQualityUber;
    params.m_endpoint_caching = false;

    dxt1_endpoint_optimizer::results results;
    results.m_pSelectors = selectors.get_ptr();

    dxt1_endpoint_optimizer optimizer;
    optimizer.compute(params, results);
    cluster.m_first_endpoint = results.m_low_color;
    cluster.m_second_endpoint = results.m_high_color;
    dxt1_block::get_block_colors4(cluster.m_color_values, cluster.m_first_endpoint, cluster.m_second_endpoint);

    color_quad_u8 color_values[4];
    color_values[0] = dxt1_block::unpack_color(results.m_low_color, true);
    color_values[3] = dxt1_block::unpack_color(results.m_high_color, true);
    for (uint c = 0; c < 3; c++) {
      color_values[1].c[c] = ((color_values[0].c[c] << 1) + color_values[3].c[c] + (results.m_alternate_rounding ? 1 : 0)) / 3;
      color_values[2].c[c] = ((color_values[3].c[c] << 1) + color_values[0].c[c] + (results.m_alternate_rounding ? 1 : 0)) / 3;
    }

    uint8 color_order[4];
    for (uint8 i = 0; i < 4; i++)
      color_order[i] = results.m_reordered ? 3 - g_dxt1_to_linear[i] : g_dxt1_to_linear[i];
    
    uint endpoint_weight = color::color_distance(m_params.m_perceptual, color_values[0], color_values[3], false) / 2000;
    float encoding_weight[8];
    for (uint i = 0; i < 8; i++)
      encoding_weight[i] = math::lerp(1.15f, 1.0f, i / 7.0f);

    crnlib::vector<uint>& blocks = cluster.m_blocks[cColorBlocks];
    for (uint i = 0; i < blocks.size(); i++) {
      uint b = blocks[i];
      uint weight = (uint)(math::clamp<uint>(endpoint_weight * m_block_weights[b], 1, 2048) * encoding_weight[m_block_encodings[b]]);
      uint32 selector = 0;
      for (uint sh = 0, p = 0; p < 16; p++, sh += 2) {
        uint error_best = cUINT32_MAX;
        uint8 s_best = 0;
        for (uint8 t = 0; t < 4; t++) {
          uint8 s = color_order[t];
          uint error = color::color_distance(m_params.m_perceptual, (color_quad_u8&)m_blocks[b][p], color_values[s], false);
          if (error < error_best) {
            s_best = s;
            error_best = error;
          }
        }
        selector |= s_best << sh;
      }
      m_block_selectors[cColorBlocks][b] = selector | (uint64)weight << 32;
    }

    dxt_endpoint_refiner refiner;
    dxt_endpoint_refiner::params refinerParams;
    dxt_endpoint_refiner::results refinerResults;
    refinerParams.m_perceptual = m_params.m_perceptual;
    refinerParams.m_pSelectors = selectors.get_ptr();
    refinerParams.m_pPixels = cluster.m_pixels.get_ptr();
    refinerParams.m_num_pixels = cluster.m_pixels.size();
    refinerParams.m_dxt1_selectors = true;
    refinerParams.m_error_to_beat = results.m_error;
    refinerParams.m_block_index = cluster_index;
    cluster.m_refined_result = refiner.refine(refinerParams, refinerResults);
    if (cluster.m_refined_result) {
      cluster.m_refined_first_endpoint = refinerResults.m_low_color;
      cluster.m_refined_second_endpoint = refinerResults.m_high_color;
    } else {
      cluster.m_refined_first_endpoint = cluster.m_first_endpoint;
      cluster.m_refined_second_endpoint = cluster.m_second_endpoint;
    }
  }
}

void dxt_hc::determine_color_endpoint_clusters_task(uint64 data, void* pData_ptr) {
  vec6F_tree_vq* vq = (vec6F_tree_vq*)pData_ptr;
  uint num_tasks = m_pTask_pool->get_num_threads() + 1;
  for (uint t = m_tiles.size() * data / num_tasks, tEnd = m_tiles.size() * (data + 1) / num_tasks; t < tEnd; t++) {
    if (m_tiles[t].pixels.size())
      m_tiles[t].cluster_indices[cColorBlocks] = vq->find_best_codebook_entry_fs(m_tiles[t].color_endpoint);
  }
}

void dxt_hc::determine_color_endpoints() {
  vec6F_tree_vq vq;
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
      m_color_clusters[m_tiles[t].cluster_indices[cColorBlocks]].m_pixels.append(m_tiles[t].pixels);
  }

  for (uint b = 0; b < m_num_blocks; b++) {
    uint cluster_index = m_tiles[m_tile_indices[b]].cluster_indices[cColorBlocks];
    m_endpoint_indices[b].component[cColorBlocks] = cluster_index;
    m_color_clusters[cluster_index].m_blocks[cColorBlocks].push_back(b);
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

    endpoint_cluster& cluster = m_alpha_clusters[cluster_index];
    if (cluster.m_pixels.empty())
      continue;

    crnlib::vector<uint8> selectors(cluster.m_pixels.size());

    dxt5_endpoint_optimizer::params params;
    params.m_pPixels = cluster.m_pixels.get_ptr();
    params.m_num_pixels = cluster.m_pixels.size();
    params.m_comp_index = 0;
    params.m_quality = cCRNDXTQualityUber;
    params.m_use_both_block_types = false;

    dxt5_endpoint_optimizer::results results;
    results.m_pSelectors = selectors.get_ptr();

    dxt5_endpoint_optimizer optimizer;
    optimizer.compute(params, results);
    cluster.m_first_endpoint = results.m_first_endpoint;
    cluster.m_second_endpoint = results.m_second_endpoint;
    dxt5_block::get_block_values(cluster.m_alpha_values, cluster.m_first_endpoint, cluster.m_second_endpoint);

    int delta = cluster.m_second_endpoint - cluster.m_first_endpoint;
    uint8 alpha_values[8];
    uint8 alpha_order[8];
    for (uint sum = cluster.m_first_endpoint * 7, i = 0; i < 8; i++, sum += delta) {
      alpha_values[i] = (uint8)(sum / 7);
      alpha_order[i] = results.m_reordered ? 7 - g_dxt5_to_linear[i] : g_dxt5_to_linear[i];
    }

    uint encoding_weight[8];
    for (uint endpoint_weight = math::clamp<uint>(delta * delta >> 3, 1, 2048), i = 0; i < 8; i++)
      encoding_weight[i] = (uint)(endpoint_weight * math::lerp(1.15f, 1.0f, i / 7.0f));

    for (uint a = 0; a < m_num_alpha_blocks; a++) {
      uint component_index = m_params.m_alpha_component_indices[a];
      crnlib::vector<uint>& blocks = cluster.m_blocks[cAlpha0Blocks + a];
      for (uint i = 0; i < blocks.size(); i++) {
        uint b = blocks[i];
        uint weight = encoding_weight[m_block_encodings[b]];
        uint64 selector = 0;
        for (uint sh = 0, p = 0; p < 16; p++, sh += 3) {
          uint error_best = cUINT32_MAX;
          uint8 s_best = 0;
          for (uint8 t = 0; t < 8; t++) {
            uint8 s = alpha_order[t];
            int delta = m_blocks[b][p][component_index] - alpha_values[s];
            uint error = delta >= 0 ? delta : -delta;
            if (error < error_best) {
              s_best = s;
              error_best = error;
            }
          }
          selector |= (uint64)s_best << sh;
        }
        m_block_selectors[cAlpha0Blocks + a][b] = selector | (uint64)weight << 48;
      }
    }

    dxt_endpoint_refiner refiner;
    dxt_endpoint_refiner::params refinerParams;
    dxt_endpoint_refiner::results refinerResults;
    refinerParams.m_perceptual = m_params.m_perceptual;
    refinerParams.m_pSelectors = selectors.get_ptr();
    refinerParams.m_pPixels = cluster.m_pixels.get_ptr();
    refinerParams.m_num_pixels = cluster.m_pixels.size();
    refinerParams.m_dxt1_selectors = false;
    refinerParams.m_error_to_beat = results.m_error;
    refinerParams.m_block_index = cluster_index;
    cluster.m_refined_result = refiner.refine(refinerParams, refinerResults);
    if (cluster.m_refined_result) {
      cluster.m_refined_first_endpoint = refinerResults.m_low_color;
      cluster.m_refined_second_endpoint = refinerResults.m_high_color;
      dxt5_block::get_block_values(cluster.m_refined_alpha_values, cluster.m_refined_first_endpoint, cluster.m_refined_second_endpoint);
    } else {
      cluster.m_refined_first_endpoint = cluster.m_first_endpoint;
      cluster.m_refined_second_endpoint = cluster.m_second_endpoint;
      memcpy(cluster.m_refined_alpha_values, cluster.m_alpha_values, sizeof(cluster.m_refined_alpha_values));
    }
  }
}

void dxt_hc::determine_alpha_endpoint_clusters_task(uint64 data, void* pData_ptr) {
  vec2F_tree_vq* vq = (vec2F_tree_vq*)pData_ptr;
  uint num_tasks = m_pTask_pool->get_num_threads() + 1;
  for (uint t = m_tiles.size() * data / num_tasks, tEnd = m_tiles.size() * (data + 1) / num_tasks; t < tEnd; t++) {
    if (m_tiles[t].pixels.size()) {
      for (uint a = 0; a < m_num_alpha_blocks; a++)
        m_tiles[t].cluster_indices[cAlpha0Blocks + a] = vq->find_best_codebook_entry_fs(m_tiles[t].alpha_endpoints[a]);
    }
  }
}

void dxt_hc::determine_alpha_endpoints() {
  vec2F_tree_vq vq;
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
        crnlib::vector<color_quad_u8>& destination = m_alpha_clusters[m_tiles[t].cluster_indices[cAlpha0Blocks + a]].m_pixels;
        for (uint p = 0; p < source.size(); p++)
          destination.push_back(color_quad_u8(source[p][component_index]));
      }
    }
  }

  for (uint b = 0; b < m_num_blocks; b++) {
    for (uint a = 0; a < m_num_alpha_blocks; a++) {
      uint cluster_index = m_tiles[m_tile_indices[b]].cluster_indices[cAlpha0Blocks + a];
      m_endpoint_indices[b].component[cAlpha0Blocks + a] = cluster_index;
      m_alpha_clusters[cluster_index].m_blocks[cAlpha0Blocks + a].push_back(b);
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
  uint errors[16][4];
  for (uint b = m_num_blocks * data / num_tasks, bEnd = m_num_blocks * (data + 1) / num_tasks; b < bEnd; b++) {
    endpoint_cluster& cluster = m_color_clusters[m_endpoint_indices[b].color];
    color_quad_u8* endpoint_colors = cluster.m_color_values;
    for (uint p = 0; p < 16; p++) {
      for (uint s = 0; s < 4; s++)
        errors[p][s] = color::color_distance(m_params.m_perceptual, m_blocks[b][p], endpoint_colors[s], false);
    }
    uint best_index = 0;
    for (uint best_error = cUINT32_MAX, s = 0; s < m_color_selectors.size(); s++) {
      uint32 selector = m_color_selectors[s];
      uint error = errors[0][selector & 3];
      error += errors[ 1][(selector >>  2) & 3];
      error += errors[ 2][(selector >>  4) & 3];
      error += errors[ 3][(selector >>  6) & 3];
      error += errors[ 4][(selector >>  8) & 3];
      error += errors[ 5][(selector >> 10) & 3];
      error += errors[ 6][(selector >> 12) & 3];
      error += errors[ 7][(selector >> 14) & 3];
      error += errors[ 8][(selector >> 16) & 3];
      error += errors[ 9][(selector >> 18) & 3];
      error += errors[10][(selector >> 20) & 3];
      error += errors[11][(selector >> 22) & 3];
      error += errors[12][(selector >> 24) & 3];
      error += errors[13][(selector >> 26) & 3];
      error += errors[14][(selector >> 28) & 3];
      error += errors[15][(selector >> 30) & 3];
      if (error < best_error) {
        best_error = error;
        best_index = s;
      }
    }
    uint (&total_errors)[16][4] = selector_details[best_index].error;
    for (uint p = 0; p < 16; p++) {
      for (uint s = 0; s < 4; s++)
        total_errors[p][s] += errors[p][s];
    }
    selector_details[best_index].used = true;
    m_selector_indices[b].color = best_index;
  }
}

void dxt_hc::create_color_selector_codebook() {
  vec16F_tree_vq selector_vq;
  vec16F v;
  for (uint b = 0; b < m_num_blocks; b++) {
    uint64 selector = m_block_selectors[cColorBlocks][b];
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
    for (uint sh = 0, j = 0; j < 16; j++, sh += 2) {
      uint8 s = g_dxt1_from_linear[(int)(v[j] * 4.0f)];
      m_color_selectors[i] |= s << sh;
    }
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
      uint best_error = errors[p][0];
      uint8 best_s = 0;
      for (uint8 s = 1; s < 4; s++) {
        uint error = errors[p][s];
        if (error < best_error) {
          best_s = s;
          best_error = error;
        }
      }
      m_color_selectors[i] |= best_s << sh;
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
  uint errors[16][8];
  for (uint b = m_num_blocks * data / num_tasks, bEnd = m_num_blocks * (data + 1) / num_tasks; b < bEnd; b++) {
    for (uint c = cAlpha0Blocks; c < cAlpha0Blocks + m_num_alpha_blocks; c++) {
      const uint alpha_pixel_comp = m_params.m_alpha_component_indices[c - cAlpha0Blocks];
      endpoint_cluster& cluster = m_alpha_clusters[m_endpoint_indices[b].component[c]];
      uint* block_values = cluster.m_alpha_values;
      for (uint p = 0; p < 16; p++) {
        for (uint s = 0; s < 8; s++) {
          int delta = m_blocks[b][p][alpha_pixel_comp] - block_values[s];
          errors[p][s] = delta * delta;
        }
      }
      uint best_index = 0;
      for (uint best_error = cUINT32_MAX, s = 0; s < m_alpha_selectors.size(); s++) {
        uint64 selector = m_alpha_selectors[s];
        uint error = errors[0][selector & 7];
        error += errors[ 1][(selector >>  3) & 7];
        error += errors[ 2][(selector >>  6) & 7];
        error += errors[ 3][(selector >>  9) & 7];
        error += errors[ 4][(selector >> 12) & 7];
        error += errors[ 5][(selector >> 15) & 7];
        error += errors[ 6][(selector >> 18) & 7];
        error += errors[ 7][(selector >> 21) & 7];
        error += errors[ 8][(selector >> 24) & 7];
        error += errors[ 9][(selector >> 27) & 7];
        error += errors[10][(selector >> 30) & 7];
        error += errors[11][(selector >> 33) & 7];
        error += errors[12][(selector >> 36) & 7];
        error += errors[13][(selector >> 39) & 7];
        error += errors[14][(selector >> 42) & 7];
        error += errors[15][(selector >> 45) & 7];
        if (error < best_error) {
          best_error = error;
          best_index = s;
        }
      }
      if (cluster.m_refined_result) {
        block_values = cluster.m_refined_alpha_values;
        for (uint p = 0; p < 16; p++) {
          for (uint s = 0; s < 8; s++) {
            int delta = m_blocks[b][p][alpha_pixel_comp] - block_values[s];
            errors[p][s] = delta * delta;
          }
        }
      }
      uint (&total_errors)[16][8] = selector_details[best_index].error;
      for (uint p = 0; p < 16; p++) {
        for (uint s = 0; s < 8; s++)
         total_errors[p][s] += errors[p][s];
      }
      selector_details[best_index].used = true;
      m_selector_indices[b].component[c] = best_index;
    }
  }
}

void dxt_hc::create_alpha_selector_codebook() {
  vec16F_tree_vq selector_vq;
  vec16F v;
  for (uint c = cAlpha0Blocks; c < cAlpha0Blocks + m_num_alpha_blocks; c++) {
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
    for (uint sh = 0, j = 0; j < 16; j++, sh += 3) {
      uint8 s = g_dxt5_from_linear[(int)(v[j] * 8.0f)];
      m_alpha_selectors[i] |= (uint64)s << sh;
    }
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
      uint best_error = errors[p][0];
      uint8 best_s = 0;
      for (uint8 s = 1; s < 8; s++) {
        uint error = errors[p][s];
        if (error < best_error) {
          best_s = s;
          best_error = error;
        }
      }
      m_alpha_selectors[i] |= (uint64)best_s << sh;
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
