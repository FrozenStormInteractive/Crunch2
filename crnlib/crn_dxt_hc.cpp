// File: crn_dxt_hc.cpp
// See Copyright Notice and license at the end of inc/crnlib.h
#include "crn_core.h"
#include "crn_dxt_hc.h"
#include "crn_image_utils.h"
#include "crn_console.h"
#include "crn_dxt_fast.h"

#define CRNLIB_USE_FAST_DXT 1
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

static color_quad_u8 g_tile_layout_colors[cNumChunkTileLayouts] =
    {
        color_quad_u8(255, 90, 32, 255),
        color_quad_u8(64, 210, 192, 255),
        color_quad_u8(128, 16, 225, 255),
        color_quad_u8(255, 192, 200, 255),

        color_quad_u8(255, 128, 200, 255),

        color_quad_u8(255, 0, 0, 255),
        color_quad_u8(0, 255, 0, 255),
        color_quad_u8(0, 0, 255, 255),
        color_quad_u8(255, 0, 255, 255)};

dxt_hc::dxt_hc()
    : m_num_chunks(0),
      m_pChunks(NULL),
      m_num_alpha_blocks(0),
      m_has_color_blocks(false),
      m_has_alpha0_blocks(false),
      m_has_alpha1_blocks(false),
      m_main_thread_id(crn_get_current_thread_id()),
      m_canceled(false),
      m_pTask_pool(NULL),
      m_prev_phase_index(-1),
      m_prev_percentage_complete(-1) {
  utils::zero_object(m_encoding_hist);
}

dxt_hc::~dxt_hc() {
}

void dxt_hc::clear() {
  m_num_chunks = 0;
  m_pChunks = NULL;

  m_num_alpha_blocks = 0;
  m_has_color_blocks = false;
  m_has_alpha0_blocks = false;
  m_has_alpha1_blocks = false;

  m_color_selectors.clear();

  m_alpha_selectors.clear();
  for (uint i = 0; i < cNumCompressedChunkVecs; i++)
    m_compressed_chunks[i].clear();

  utils::zero_object(m_encoding_hist);

  m_total_tiles = 0;

  m_color_clusters.clear();
  m_alpha_clusters.clear();
  m_color_selectors.clear();
  m_alpha_selectors.clear();

  m_chunk_blocks_using_color_selectors.clear();
  m_chunk_blocks_using_alpha_selectors.clear();

  m_color_endpoints.clear();
  m_alpha_endpoints.clear();

  m_canceled = false;

  m_prev_phase_index = -1;
  m_prev_percentage_complete = -1;

  m_chunk_details.clear();
  m_blocks.clear();
  for (uint c = 0; c < 3; c++)
    m_block_selectors[c].clear();
  m_endpoint_indices.clear();

}

bool dxt_hc::initialize_blocks(const params& p) {
  m_chunk_details.resize(m_num_chunks);
  m_blocks.resize(m_num_chunks << 2);
  for (uint c = 0; c < 3; c++)
    m_block_selectors[c].resize(m_blocks.size());
  m_endpoint_indices.resize(m_blocks.size());

  for (uint level = 0; level < p.m_num_levels; level++) {
    uint first_chunk = p.m_levels[level].m_first_chunk;
    uint end_chunk = p.m_levels[level].m_first_chunk + p.m_levels[level].m_num_chunks;
    uint chunk_width = p.m_levels[level].m_chunk_width;
    uint block_width = chunk_width << 1;
    for (uint b = first_chunk << 2, cy = 0, chunk_base = first_chunk; chunk_base < end_chunk; chunk_base += chunk_width, cy++) {
      for (uint by = 0; by < 2; by++) {
        for (uint cx = 0; cx < chunk_width; cx++) {
          for (uint bx = 0; bx < 2; bx++, b++) {
            const pixel_chunk& chunk = m_pChunks[chunk_base + cx];
            m_chunk_details[chunk_base + cx].block_index[by][bx] = b;
            for (uint t = 0, y = 0; y < 4; y++)  {
              for (uint x = 0; x < 4; x++, t++)
                m_blocks[b].push_back(chunk(bx << 2 | x, by << 2 | y));
            }
          }
        }
      }
    }
  }
  return true;
}

bool dxt_hc::compress(const params& p, uint num_chunks, const pixel_chunk* pChunks, task_pool& task_pool) {
  m_pTask_pool = &task_pool;
  m_main_thread_id = crn_get_current_thread_id();

  bool result = compress_internal(p, num_chunks, pChunks);

  m_pTask_pool = NULL;

  return result;
}

bool dxt_hc::compress_internal(const params& p, uint num_chunks, const pixel_chunk* pChunks) {
  if ((!num_chunks) || (!pChunks))
    return false;
  if ((m_params.m_format == cDXT1A) || (m_params.m_format == cDXT3))
    return false;

  clear();

  m_params = p;

  m_num_chunks = num_chunks;
  m_pChunks = pChunks;

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

  initialize_blocks(p);
  determine_compressed_chunks();

  if (m_has_color_blocks) {
    if (!determine_color_endpoint_clusters())
      return false;
    if (!determine_color_endpoint_codebook())
      return false;
  }

  if (m_num_alpha_blocks) {
    if (!determine_alpha_endpoint_clusters())
      return false;
    if (!determine_alpha_endpoint_codebook())
      return false;
  }

  if (m_has_color_blocks) {
    if (!create_selector_codebook(false))
      return false;
  }

  if (m_num_alpha_blocks) {
    if (!create_selector_codebook(true))
      return false;
  }

  if (m_has_color_blocks) {
    if (!refine_quantized_color_selectors())
      return false;

    if (!refine_quantized_color_endpoints())
      return false;
  }

  if (m_num_alpha_blocks) {
    if (!refine_quantized_alpha_endpoints())
      return false;

    if (!refine_quantized_alpha_selectors())
      return false;
  }

  if (!create_block_encodings(p))
    return false;

  return true;
}

void dxt_hc::compress_dxt1_block(
    dxt1_endpoint_optimizer::results& results,
    uint chunk_index, const image_u8& chunk, uint x_ofs, uint y_ofs, uint width, uint height,
    uint8* pColor_Selectors) {
  chunk_index;

  color_quad_u8 pixels[cChunkPixelWidth * cChunkPixelHeight];

  for (uint y = 0; y < height; y++)
    for (uint x = 0; x < width; x++)
      pixels[x + y * width] = chunk(x_ofs + x, y_ofs + y);

//double s = image_utils::compute_std_dev(width * height, pixels, 0, 3);

#if CRNLIB_USE_FAST_DXT
  uint low16, high16;
  dxt_fast::compress_color_block(width * height, pixels, low16, high16, pColor_Selectors);
  results.m_low_color = static_cast<uint16>(low16);
  results.m_high_color = static_cast<uint16>(high16);
  results.m_alpha_block = false;
  results.m_error = INT_MAX;
  results.m_pSelectors = pColor_Selectors;
#else
  dxt1_endpoint_optimizer optimizer;

  dxt1_endpoint_optimizer::params params;
  params.m_block_index = chunk_index;
  params.m_pPixels = pixels;
  params.m_num_pixels = width * height;
  params.m_pixels_have_alpha = false;
  params.m_use_alpha_blocks = false;
  params.m_perceptual = m_params.m_perceptual;
  params.m_highest_quality = false;  //false;
  params.m_endpoint_caching = false;

  results.m_pSelectors = pColor_Selectors;

  optimizer.compute(params, results);
#endif
}

void dxt_hc::compress_dxt5_block(
    dxt5_endpoint_optimizer::results& results,
    uint chunk_index, const image_u8& chunk, uint x_ofs, uint y_ofs, uint width, uint height, uint component_index,
    uint8* pAlpha_selectors) {
  chunk_index;

  color_quad_u8 pixels[cChunkPixelWidth * cChunkPixelHeight];

  for (uint y = 0; y < height; y++)
    for (uint x = 0; x < width; x++)
      pixels[x + y * width] = chunk(x_ofs + x, y_ofs + y);

#if 0  //CRNLIB_USE_FAST_DXT
      uint low, high;
      dxt_fast::compress_alpha_block(width * height, pixels, low, high, pAlpha_selectors, component_index);
      results.m_pSelectors = pAlpha_selectors;
      results.m_error = INT_MAX;
      results.m_first_endpoint = static_cast<uint8>(low);
      results.m_second_endpoint = static_cast<uint8>(high);
      results.m_block_type = 0;
#else
  dxt5_endpoint_optimizer optimizer;
  dxt5_endpoint_optimizer::params params;
  params.m_block_index = chunk_index;
  params.m_pPixels = pixels;
  params.m_num_pixels = width * height;
  params.m_comp_index = component_index;
  params.m_use_both_block_types = false;
  params.m_quality = cCRNDXTQualityNormal;

  results.m_pSelectors = pAlpha_selectors;

  optimizer.compute(params, results);
#endif
}

void dxt_hc::determine_compressed_chunks_task(uint64 data, void* pData_ptr) {
  pData_ptr;
  const uint thread_index = static_cast<uint>(data);

  image_u8 orig_chunk;
  image_u8 decomp_chunk[cNumChunkEncodings];

  orig_chunk.resize(cChunkPixelWidth, cChunkPixelHeight);
  for (uint i = 0; i < cNumChunkEncodings; i++)
    decomp_chunk[i].resize(cChunkPixelWidth, cChunkPixelHeight);

  image_utils::error_metrics color_error_metrics[cNumChunkEncodings];
  dxt1_endpoint_optimizer::results color_optimizer_results[cNumChunkTileLayouts];
  uint8 layout_color_selectors[cNumChunkTileLayouts][cChunkPixelWidth * cChunkPixelHeight];

  image_utils::error_metrics alpha_error_metrics[2][cNumChunkEncodings];
  dxt5_endpoint_optimizer::results alpha_optimizer_results[2][cNumChunkTileLayouts];
  uint8 layout_alpha_selectors[2][cNumChunkTileLayouts][cChunkPixelWidth * cChunkPixelHeight];

  uint first_layout = 0;
  uint last_layout = cNumChunkTileLayouts;

  uint first_encoding = 0;
  uint last_encoding = cNumChunkEncodings;

  if (!m_params.m_hierarchical) {
    first_layout = cFirst4x4ChunkTileLayout;
    first_encoding = cNumChunkEncodings - 1;
  }

  for (uint chunk_index = 0; chunk_index < m_num_chunks; chunk_index++) {
    if (m_canceled)
      return;

    if ((crn_get_current_thread_id() == m_main_thread_id) && ((chunk_index & 511) == 0)) {
      if (!update_progress(0, chunk_index, m_num_chunks))
        return;
    }

    if (m_pTask_pool->get_num_threads()) {
      if ((chunk_index % (m_pTask_pool->get_num_threads() + 1)) != thread_index)
        continue;
    }

    uint level_index = 0;
    for (uint i = 0; i < m_params.m_num_levels; i++) {
      if ((chunk_index >= m_params.m_levels[i].m_first_chunk) && (chunk_index < m_params.m_levels[i].m_first_chunk + m_params.m_levels[i].m_num_chunks)) {
        level_index = i;
        break;
      }
    }

    for (uint cy = 0; cy < cChunkPixelHeight; cy++)
      for (uint cx = 0; cx < cChunkPixelWidth; cx++)
        orig_chunk(cx, cy) = m_pChunks[chunk_index](cx, cy);

    if (m_has_color_blocks) {
      for (uint l = first_layout; l < last_layout; l++) {
        utils::zero_object(layout_color_selectors[l]);

        compress_dxt1_block(
            color_optimizer_results[l], chunk_index,
            orig_chunk,
            g_chunk_tile_layouts[l].m_x_ofs, g_chunk_tile_layouts[l].m_y_ofs,
            g_chunk_tile_layouts[l].m_width, g_chunk_tile_layouts[l].m_height,
            layout_color_selectors[l]);
      }
    }

    float alpha_layout_std_dev[2][cNumChunkTileLayouts];
    utils::zero_object(alpha_layout_std_dev);

    for (uint a = 0; a < m_num_alpha_blocks; a++) {
      for (uint l = first_layout; l < last_layout; l++) {
        utils::zero_object(layout_alpha_selectors[a][l]);

        compress_dxt5_block(
            alpha_optimizer_results[a][l], chunk_index,
            orig_chunk,
            g_chunk_tile_layouts[l].m_x_ofs, g_chunk_tile_layouts[l].m_y_ofs,
            g_chunk_tile_layouts[l].m_width, g_chunk_tile_layouts[l].m_height,
            m_params.m_alpha_component_indices[a],
            layout_alpha_selectors[a][l]);

        for (uint a = 0; a < m_num_alpha_blocks; a++) {
          float mean = 0.0f;
          float variance = 0.0f;

          for (uint cy = 0; cy < g_chunk_tile_layouts[l].m_height; cy++) {
            for (uint cx = 0; cx < g_chunk_tile_layouts[l].m_width; cx++) {
              uint s = orig_chunk(cx + g_chunk_tile_layouts[l].m_x_ofs, cy + g_chunk_tile_layouts[l].m_y_ofs)[m_params.m_alpha_component_indices[a]];

              mean += s;
              variance += s * s;
            }  // cx
          }    //cy

          float scale = 1.0f / (g_chunk_tile_layouts[l].m_width * g_chunk_tile_layouts[l].m_height);

          mean *= scale;
          variance *= scale;

          variance -= mean * mean;

          alpha_layout_std_dev[a][l] = sqrt(variance);

        }  //a
      }
    }

    for (uint e = first_encoding; e < last_encoding; e++) {
      for (uint t = 0; t < g_chunk_encodings[e].m_num_tiles; t++) {
        const uint layout_index = g_chunk_encodings[e].m_tiles[t].m_layout_index;
        CRNLIB_ASSERT((layout_index >= first_layout) && (layout_index < last_layout));

        if (m_has_color_blocks) {
          const dxt1_endpoint_optimizer::results& color_results = color_optimizer_results[layout_index];
          const uint8* pColor_selectors = layout_color_selectors[layout_index];

          color_quad_u8 block_colors[cDXT1SelectorValues];
          CRNLIB_ASSERT(color_results.m_low_color >= color_results.m_high_color);
          // it's okay if color_results.m_low_color == color_results.m_high_color, because in this case only selector 0 should be used
          dxt1_block::get_block_colors4(block_colors, color_results.m_low_color, color_results.m_high_color);

          for (uint cy = 0; cy < g_chunk_encodings[e].m_tiles[t].m_height; cy++) {
            for (uint cx = 0; cx < g_chunk_encodings[e].m_tiles[t].m_width; cx++) {
              uint s = pColor_selectors[cx + cy * g_chunk_encodings[e].m_tiles[t].m_width];
              CRNLIB_ASSERT(s < cDXT1SelectorValues);

              decomp_chunk[e](cx + g_chunk_encodings[e].m_tiles[t].m_x_ofs, cy + g_chunk_encodings[e].m_tiles[t].m_y_ofs) = block_colors[s];
            }
          }
        }

        for (uint a = 0; a < m_num_alpha_blocks; a++) {
          const dxt5_endpoint_optimizer::results& alpha_results = alpha_optimizer_results[a][layout_index];
          const uint8* pAlpha_selectors = layout_alpha_selectors[a][layout_index];

          uint block_values[cDXT5SelectorValues];
          CRNLIB_ASSERT(alpha_results.m_first_endpoint >= alpha_results.m_second_endpoint);
          dxt5_block::get_block_values8(block_values, alpha_results.m_first_endpoint, alpha_results.m_second_endpoint);

          for (uint cy = 0; cy < g_chunk_encodings[e].m_tiles[t].m_height; cy++) {
            for (uint cx = 0; cx < g_chunk_encodings[e].m_tiles[t].m_width; cx++) {
              uint s = pAlpha_selectors[cx + cy * g_chunk_encodings[e].m_tiles[t].m_width];
              CRNLIB_ASSERT(s < cDXT5SelectorValues);

              decomp_chunk[e](cx + g_chunk_encodings[e].m_tiles[t].m_x_ofs, cy + g_chunk_encodings[e].m_tiles[t].m_y_ofs)[m_params.m_alpha_component_indices[a]] =
                  static_cast<uint8>(block_values[s]);
            }
          }
        }
      }  // t

      if (m_params.m_hierarchical) {
        if (m_has_color_blocks)
          color_error_metrics[e].compute(decomp_chunk[e], orig_chunk, 0, 3);

        for (uint a = 0; a < m_num_alpha_blocks; a++)
          alpha_error_metrics[a][e].compute(decomp_chunk[e], orig_chunk, m_params.m_alpha_component_indices[a], 1);
      }
    }  // e

    uint best_encoding = cNumChunkEncodings - 1;

    if (m_params.m_hierarchical) {
      float quality[cNumChunkEncodings];
      utils::zero_object(quality);

      float best_quality = 0.0f;

      best_encoding = 0;

      for (uint e = 0; e < cNumChunkEncodings; e++) {
        if (m_has_color_blocks) {
          float adaptive_tile_color_psnr_derating = m_params.m_adaptive_tile_color_psnr_derating;
          if ((level_index) && (adaptive_tile_color_psnr_derating > .25f)) {
            //adaptive_tile_color_psnr_derating = math::lerp(adaptive_tile_color_psnr_derating * .5f, .3f, (level_index - 1) / math::maximum(1.0f, float(m_params.m_num_levels - 2)));
            adaptive_tile_color_psnr_derating = math::maximum(.25f, adaptive_tile_color_psnr_derating / powf(3.0f, static_cast<float>(level_index)));
          }

          float color_derating = math::lerp(0.0f, adaptive_tile_color_psnr_derating, (g_chunk_encodings[e].m_num_tiles - 1) / 3.0f);
          quality[e] = (float)math::maximum<double>(color_error_metrics[e].mPeakSNR - color_derating, 0.0f);
        }

        if (m_num_alpha_blocks) {
          quality[e] *= m_params.m_adaptive_tile_color_alpha_weighting_ratio;
          float alpha_derating = math::lerp(0.0f, m_params.m_adaptive_tile_alpha_psnr_derating, (g_chunk_encodings[e].m_num_tiles - 1) / 3.0f);

          float max_std_dev = 0.0f;

          for (uint a = 0; a < m_num_alpha_blocks; a++) {
            quality[e] += (float)math::maximum<double>(alpha_error_metrics[a][e].mPeakSNR - alpha_derating, 0.0f);

            for (uint t = 0; t < g_chunk_encodings[e].m_num_tiles; t++) {
              float std_dev = alpha_layout_std_dev[a][g_chunk_encodings[e].m_tiles[t].m_layout_index];
              max_std_dev = math::maximum(max_std_dev, std_dev);
            }
          }

#if 0
// rg [4/28/09] - disabling this because it's fucking up dxt5_xgbr normal maps
                  const float l = 6.0f;
                  const float k = .5f;

                  if (max_std_dev > l)
                  {
                     float s = max_std_dev - l;
                     quality[e] -= (k * s);
                  }
#endif
        }

        if (quality[e] > best_quality) {
          best_quality = quality[e];
          best_encoding = e;
        }
      }
    }

    atomic_increment32(&m_encoding_hist[best_encoding]);

    atomic_exchange_add32(&m_total_tiles, g_chunk_encodings[best_encoding].m_num_tiles);

    for (uint q = 0; q < cNumCompressedChunkVecs; q++) {
      if (q == cColorChunks) {
        if (!m_has_color_blocks)
          continue;
      } else if (q > m_num_alpha_blocks)
        continue;

      compressed_chunk& output = m_compressed_chunks[q][chunk_index];

      output.m_encoding_index = static_cast<uint8>(best_encoding);
      output.m_num_tiles = static_cast<uint8>(g_chunk_encodings[best_encoding].m_num_tiles);

      for (uint t = 0; t < g_chunk_encodings[best_encoding].m_num_tiles; t++) {
        const uint layout_index = g_chunk_encodings[best_encoding].m_tiles[t].m_layout_index;

        output.m_tiles[t].m_layout_index = static_cast<uint8>(layout_index);
        output.m_tiles[t].m_pixel_width = static_cast<uint8>(g_chunk_encodings[best_encoding].m_tiles[t].m_width);
        output.m_tiles[t].m_pixel_height = static_cast<uint8>(g_chunk_encodings[best_encoding].m_tiles[t].m_height);

        if (q == cColorChunks) {
          const dxt1_endpoint_optimizer::results& color_results = color_optimizer_results[layout_index];
          const uint8* pColor_selectors = layout_color_selectors[layout_index];

          output.m_tiles[t].m_endpoint_cluster_index = 0;
          output.m_tiles[t].m_first_endpoint = color_results.m_low_color;
          output.m_tiles[t].m_second_endpoint = color_results.m_high_color;

        } else {
          const uint a = q - cAlpha0Chunks;

          const dxt5_endpoint_optimizer::results& alpha_results = alpha_optimizer_results[a][layout_index];
          const uint8* pAlpha_selectors = layout_alpha_selectors[a][layout_index];

          output.m_tiles[t].m_endpoint_cluster_index = 0;
          output.m_tiles[t].m_first_endpoint = alpha_results.m_first_endpoint;
          output.m_tiles[t].m_second_endpoint = alpha_results.m_second_endpoint;

        }
      }  // t
    }    // q
  }  // chunk_index
}

bool dxt_hc::determine_compressed_chunks() {
  utils::zero_object(m_encoding_hist);

  for (uint i = 0; i < cNumCompressedChunkVecs; i++)
    m_compressed_chunks[i].clear();

  if (m_has_color_blocks)
    m_compressed_chunks[cColorChunks].resize(m_num_chunks);

  for (uint a = 0; a < m_num_alpha_blocks; a++)
    m_compressed_chunks[cAlpha0Chunks + a].resize(m_num_chunks);

  m_total_tiles = 0;

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::determine_compressed_chunks_task, i);

  m_pTask_pool->join();
  if (m_canceled)
    return false;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging) {
    console::info("Total Pixels: %u, Chunks: %u, Blocks: %u, Adapted Tiles: %u", m_num_chunks * cChunkPixelWidth * cChunkPixelHeight, m_num_chunks, m_num_chunks * cChunkBlockWidth * cChunkBlockHeight, m_total_tiles);

    console::info("Chunk encoding type symbol_histogram: ");
    for (uint e = 0; e < cNumChunkEncodings; e++)
      console::info("%u ", m_encoding_hist[e]);

    console::info("Blocks per chunk encoding type: ");
    for (uint e = 0; e < cNumChunkEncodings; e++)
      console::info("%u ", m_encoding_hist[e] * cChunkBlockWidth * cChunkBlockHeight);
  }
#endif

  return true;
}

void dxt_hc::assign_color_endpoint_clusters_task(uint64 data, void* pData_ptr) {
  const uint thread_index = (uint)data;
  assign_color_endpoint_clusters_state& state = *static_cast<assign_color_endpoint_clusters_state*>(pData_ptr);

  for (uint chunk_index = 0; chunk_index < m_num_chunks; chunk_index++) {
    if (m_canceled)
      return;

    if ((crn_get_current_thread_id() == m_main_thread_id) && ((chunk_index & 63) == 0)) {
      if (!update_progress(2, chunk_index, m_num_chunks))
        return;
    }

    if (m_pTask_pool->get_num_threads()) {
      if ((chunk_index % (m_pTask_pool->get_num_threads() + 1)) != thread_index)
        continue;
    }

    compressed_chunk& chunk = m_compressed_chunks[cColorChunks][chunk_index];

    for (uint tile_index = 0; tile_index < chunk.m_num_tiles; tile_index++) {
      uint cluster_index = state.m_vq.find_best_codebook_entry_fs(state.m_training_vecs[chunk_index][tile_index]);

      chunk.m_endpoint_cluster_index[tile_index] = static_cast<uint16>(cluster_index);
    }
  }
}

bool dxt_hc::determine_color_endpoint_clusters() {
  if (!m_has_color_blocks)
    return true;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Generating color training vectors");
#endif

  const float r_scale = .5f;
  const float b_scale = .25f;

  vec6F_tree_vq vq;

  crnlib::vector<crnlib::vector<vec6F> > training_vecs;

  training_vecs.resize(m_num_chunks);

  for (uint chunk_index = 0; chunk_index < m_num_chunks; chunk_index++) {
    if ((chunk_index & 255) == 0) {
      if (!update_progress(1, chunk_index, m_num_chunks))
        return false;
    }

    const compressed_chunk& chunk = m_compressed_chunks[cColorChunks][chunk_index];

    training_vecs[chunk_index].resize(chunk.m_num_tiles);

    for (uint tile_index = 0; tile_index < chunk.m_num_tiles; tile_index++) {
      const compressed_tile& tile = chunk.m_tiles[tile_index];

      const chunk_tile_desc& layout = g_chunk_tile_layouts[tile.m_layout_index];

      tree_clusterizer<vec3F> palettizer;
      for (uint y = 0; y < layout.m_height; y++) {
        for (uint x = 0; x < layout.m_width; x++) {
          const color_quad_u8& c = m_pChunks[chunk_index](layout.m_x_ofs + x, layout.m_y_ofs + y);

          vec3F v;
          if (m_params.m_perceptual) {
            v.set(c[0] * 1.0f / 255.0f, c[1] * 1.0f / 255.0f, c[2] * 1.0f / 255.0f);
            v[0] *= r_scale;
            v[2] *= b_scale;
          } else {
            v.set(c[0] * 1.0f / 255.0f, c[1] * 1.0f / 255.0f, c[2] * 1.0f / 255.0f);
          }

          palettizer.add_training_vec(v, 1);
        }
      }

      palettizer.generate_codebook(2);

      uint tile_weight = tile.m_pixel_width * tile.m_pixel_height;
      tile_weight = static_cast<uint>(tile_weight * m_pChunks[chunk_index].m_weight);

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

      vq.add_training_vec(vv, tile_weight);

      training_vecs[chunk_index][tile_index] = vv;
    }
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Begin color cluster analysis");
  timer t;
  t.start();
#endif

  uint codebook_size = math::minimum<uint>(m_total_tiles, m_params.m_color_endpoint_codebook_size);
  vq.generate_codebook(codebook_size);

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging) {
    double total_time = t.get_elapsed_secs();
    console::info("Codebook gen time: %3.3fs, Total color clusters: %u", total_time, vq.get_codebook_size());
  }
#endif

  m_color_clusters.resize(vq.get_codebook_size());

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Begin color cluster assignment");
#endif

  assign_color_endpoint_clusters_state state(vq, training_vecs);

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::assign_color_endpoint_clusters_task, i, &state);

  m_pTask_pool->join();
  if (m_canceled)
    return false;

  for (uint i = 0; i < m_num_chunks; i++) {
    int chunk_index = m_pChunks[i].m_legacy_index;
    compressed_chunk& chunk = m_compressed_chunks[cColorChunks][chunk_index];
    for (uint tile_index = 0; tile_index < chunk.m_num_tiles; tile_index++) {
      uint cluster_index = chunk.m_endpoint_cluster_index[tile_index];
      m_color_clusters[cluster_index].m_tiles.push_back(std::make_pair(chunk_index, tile_index));

      const compressed_tile& tile = chunk.m_tiles[tile_index];
      const chunk_tile_desc& layout = g_chunk_tile_layouts[tile.m_layout_index];
      for (uint y = 0; y < layout.m_height; y++)
        for (uint x = 0; x < layout.m_width; x++)
          m_color_clusters[cluster_index].m_pixels.push_back(m_pChunks[chunk_index](layout.m_x_ofs + x, layout.m_y_ofs + y));
    }
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Completed color cluster assignment");
#endif

  return true;
}

void dxt_hc::determine_alpha_endpoint_clusters_task(uint64 data, void* pData_ptr) {
  const uint thread_index = static_cast<uint>(data);
  const determine_alpha_endpoint_clusters_state& state = *static_cast<determine_alpha_endpoint_clusters_state*>(pData_ptr);

  for (uint a = 0; a < m_num_alpha_blocks; a++) {
    for (uint chunk_index = 0; chunk_index < m_num_chunks; chunk_index++) {
      if (m_canceled)
        return;

      if ((crn_get_current_thread_id() == m_main_thread_id) && ((chunk_index & 63) == 0)) {
        if (!update_progress(7, m_num_chunks * a + chunk_index, m_num_chunks * m_num_alpha_blocks))
          return;
      }

      if (m_pTask_pool->get_num_threads()) {
        if ((chunk_index % (m_pTask_pool->get_num_threads() + 1)) != thread_index)
          continue;
      }

      compressed_chunk& chunk = m_compressed_chunks[cAlpha0Chunks + a][chunk_index];

      for (uint tile_index = 0; tile_index < chunk.m_num_tiles; tile_index++) {
        uint cluster_index = state.m_vq.find_best_codebook_entry_fs(state.m_training_vecs[a][chunk_index][tile_index]);

        chunk.m_endpoint_cluster_index[tile_index] = static_cast<uint16>(cluster_index);
      }
    }
  }
}

bool dxt_hc::determine_alpha_endpoint_clusters() {
  if (!m_num_alpha_blocks)
    return true;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Generating alpha training vectors");
#endif

  determine_alpha_endpoint_clusters_state state;

  for (uint a = 0; a < m_num_alpha_blocks; a++) {
    state.m_training_vecs[a].resize(m_num_chunks);

    for (uint chunk_index = 0; chunk_index < m_num_chunks; chunk_index++) {
      if ((chunk_index & 63) == 0) {
        if (!update_progress(6, m_num_chunks * a + chunk_index, m_num_chunks * m_num_alpha_blocks))
          return false;
      }

      const compressed_chunk& chunk = m_compressed_chunks[cAlpha0Chunks + a][chunk_index];

      state.m_training_vecs[a][chunk_index].resize(chunk.m_num_tiles);

      for (uint tile_index = 0; tile_index < chunk.m_num_tiles; tile_index++) {
        const compressed_tile& tile = chunk.m_tiles[tile_index];

        const chunk_tile_desc& layout = g_chunk_tile_layouts[tile.m_layout_index];

        tree_clusterizer<vec1F> palettizer;

        for (uint y = 0; y < layout.m_height; y++) {
          for (uint x = 0; x < layout.m_width; x++) {
            uint c = m_pChunks[chunk_index](layout.m_x_ofs + x, layout.m_y_ofs + y)[m_params.m_alpha_component_indices[a]];

            vec1F v(c * 1.0f / 255.0f);

            palettizer.add_training_vec(v, 1);
          }
        }
        palettizer.generate_codebook(2);

        const uint tile_weight = tile.m_pixel_width * tile.m_pixel_height;

        vec1F v[2];
        utils::zero_object(v);

        for (uint i = 0; i < palettizer.get_codebook_size(); i++)
          v[i] = palettizer.get_codebook_entry(i);

        if (palettizer.get_codebook_size() == 1)
          v[1] = v[0];
        if (v[0] > v[1])
          utils::swap(v[0], v[1]);

        vec2F vv(v[0][0], v[1][0]);

        state.m_vq.add_training_vec(vv, tile_weight);

        state.m_training_vecs[a][chunk_index][tile_index] = vv;

      }  // tile_index
    }    // chunk_index
  }      // a

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Begin alpha cluster analysis");
  timer t;
  t.start();
#endif

  uint codebook_size = math::minimum<uint>(m_total_tiles, m_params.m_alpha_endpoint_codebook_size);
  state.m_vq.generate_codebook(codebook_size);

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging) {
    double total_time = t.get_elapsed_secs();
    console::info("Codebook gen time: %3.3fs, Total alpha clusters: %u", total_time, state.m_vq.get_codebook_size());
  }
#endif

  m_alpha_clusters.resize(state.m_vq.get_codebook_size());

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Begin alpha cluster assignment");
#endif

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::determine_alpha_endpoint_clusters_task, i, &state);

  m_pTask_pool->join();
  if (m_canceled)
    return false;

  for (uint a = 0; a < m_num_alpha_blocks; a++) {
    uint component_index = m_params.m_alpha_component_indices[a];
    for (uint i = 0; i < m_num_chunks; i++) {
      int chunk_index = m_pChunks[i].m_legacy_index;
      compressed_chunk& chunk = m_compressed_chunks[cAlpha0Chunks + a][chunk_index];
      for (uint tile_index = 0; tile_index < chunk.m_num_tiles; tile_index++) {
        const uint cluster_index = chunk.m_endpoint_cluster_index[tile_index];
        m_alpha_clusters[cluster_index].m_tiles.push_back(std::make_pair(chunk_index, tile_index | (a << 16)));

        const compressed_tile& tile = chunk.m_tiles[tile_index];
        const chunk_tile_desc& layout = g_chunk_tile_layouts[tile.m_layout_index];
        for (uint y = 0; y < layout.m_height; y++)
          for (uint x = 0; x < layout.m_width; x++)
            m_alpha_clusters[cluster_index].m_pixels.push_back(color_quad_u8(m_pChunks[chunk_index](layout.m_x_ofs + x, layout.m_y_ofs + y)[component_index]));
      }
    }
  }

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Completed alpha cluster assignment");
#endif

  return true;
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

    tile_cluster& cluster = m_color_clusters[cluster_index];
    if (cluster.m_pixels.empty())
      continue;

    cluster.m_selectors.resize(cluster.m_pixels.size());

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
    results.m_pSelectors = cluster.m_selectors.get_ptr();

    dxt1_endpoint_optimizer optimizer;
    optimizer.compute(params, results);

    cluster.m_first_endpoint = results.m_low_color;
    cluster.m_second_endpoint = results.m_high_color;
    cluster.m_error = results.m_error;

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

    for (uint t = 0; t < cluster.m_tiles.size(); t++) {
      const uint chunk_index = cluster.m_tiles[t].first;
      const uint tile_index = cluster.m_tiles[t].second;
      compressed_chunk& chunk = m_compressed_chunks[cColorChunks][chunk_index];
      uint8 encoding_index = chunk.m_encoding_index;
      uint weight = (uint)(math::clamp<uint>(endpoint_weight * m_pChunks[chunk_index].m_weight, 1, 2048) * encoding_weight[encoding_index]);
      for (uint by = 0; by < 2; by++) {
        for (uint bx = 0; bx < 2; bx++) {
          if (g_tile_map[encoding_index][by][bx] == tile_index) {
            uint b = m_chunk_details[chunk_index].block_index[by][bx];
            uint64 selector = 0;
            for (uint sh = 0, p = 0; p < 16; p++, sh += 3) {
              uint8 s_best;
              for (uint32 error_best = UINT_MAX, t = 0; t < 4; t++) {
                uint8 s = color_order[t];
                uint32 error = color::color_distance(m_params.m_perceptual, (color_quad_u8&)m_blocks[b][p], color_values[s], false);
                if (error < error_best) {
                  s_best = s;
                  error_best = error;
                }
              }
              selector |= (uint64)s_best << sh;
            }
            m_block_selectors[cColorChunks][b] = selector | (uint64)weight << 48;
            m_endpoint_indices[b].component[0] = cluster_index;
          }
        }
      }
    }

    dxt_endpoint_refiner refiner;
    dxt_endpoint_refiner::params refinerParams;
    dxt_endpoint_refiner::results refinerResults;
    refinerParams.m_perceptual = m_params.m_perceptual;
    refinerParams.m_pSelectors = cluster.m_selectors.get_ptr();
    refinerParams.m_pPixels = cluster.m_pixels.get_ptr();
    refinerParams.m_num_pixels = cluster.m_pixels.size();
    refinerParams.m_dxt1_selectors = true;
    refinerParams.m_error_to_beat = cluster.m_error;
    refinerParams.m_block_index = cluster_index;
    cluster.m_refined.result = refiner.refine(refinerParams, refinerResults);
    cluster.m_refined.first_endpoint = refinerResults.m_low_color;
    cluster.m_refined.second_endpoint = refinerResults.m_high_color;
    cluster.m_refined.error = refinerResults.m_error;

    for (uint t = 0; t < cluster.m_tiles.size(); t++) {
      const uint chunk_index = cluster.m_tiles[t].first;
      const uint tile_index = cluster.m_tiles[t].second;

      CRNLIB_ASSERT(chunk_index < m_num_chunks);

      compressed_chunk& chunk = m_compressed_chunks[cColorChunks][chunk_index];

      CRNLIB_ASSERT(tile_index < chunk.m_num_tiles);

      CRNLIB_ASSERT(chunk.m_endpoint_cluster_index[tile_index] == cluster_index);

      const compressed_tile& tile = chunk.m_tiles[tile_index];

      const chunk_tile_desc& layout = g_chunk_tile_layouts[tile.m_layout_index];
      layout;

      compressed_tile& quantized_tile = chunk.m_quantized_tiles[tile_index];

      const uint total_pixels = tile.m_pixel_width * tile.m_pixel_height;

      quantized_tile.m_endpoint_cluster_index = cluster_index;
      quantized_tile.m_first_endpoint = results.m_low_color;
      quantized_tile.m_second_endpoint = results.m_high_color;
      quantized_tile.m_pixel_width = tile.m_pixel_width;
      quantized_tile.m_pixel_height = tile.m_pixel_height;
      quantized_tile.m_layout_index = tile.m_layout_index;
    }
  }
}

bool dxt_hc::determine_color_endpoint_codebook() {
  if (!m_has_color_blocks)
    return true;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Computing optimal color cluster endpoints");
#endif

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::determine_color_endpoint_codebook_task, i, NULL);

  m_pTask_pool->join();

  return !m_canceled;
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

    tile_cluster& cluster = m_alpha_clusters[cluster_index];
    if (cluster.m_pixels.empty())
      continue;

    cluster.m_selectors.resize(cluster.m_pixels.size());

    dxt5_endpoint_optimizer::params params;
    params.m_block_index = cluster_index;
    params.m_pPixels = cluster.m_pixels.get_ptr();
    params.m_num_pixels = cluster.m_pixels.size();
    params.m_comp_index = 0;
    params.m_quality = cCRNDXTQualityUber;
    params.m_use_both_block_types = false;

    dxt5_endpoint_optimizer::results results;
    results.m_pSelectors = cluster.m_selectors.get_ptr();

    dxt5_endpoint_optimizer optimizer;
    optimizer.compute(params, results);

    cluster.m_first_endpoint = results.m_first_endpoint;
    cluster.m_second_endpoint = results.m_second_endpoint;
    cluster.m_error = results.m_error;

    int delta = cluster.m_second_endpoint - cluster.m_first_endpoint;
    uint8 alpha_values[8];
    uint8 alpha_order[8];
    for (uint sum = cluster.m_first_endpoint * 7, i = 0; i < 8; i++, sum += delta) {
      alpha_values[i] = (uint8)(sum / 7);
      alpha_order[i] = results.m_reordered ? 7 - g_dxt5_to_linear[i] : g_dxt5_to_linear[i];
    }
    uint64 encoding_weight[8];
    for (uint endpoint_weight = math::clamp<uint>(delta * delta >> 3, 1, 2048), i = 0; i < 8; i++)
      encoding_weight[i] = (uint)(endpoint_weight * math::lerp(1.15f, 1.0f, i / 7.0f));

    for (uint tile_iter = 0; tile_iter < cluster.m_tiles.size(); tile_iter++) {
      const uint chunk_index = cluster.m_tiles[tile_iter].first;
      const uint tile_index = cluster.m_tiles[tile_iter].second & 0xFFFFU;
      const uint alpha_index = cluster.m_tiles[tile_iter].second >> 16U;
      compressed_chunk& chunk = m_compressed_chunks[cAlpha0Chunks + alpha_index][chunk_index];
      uint component_index = m_params.m_alpha_component_indices[alpha_index];
      uint8 encoding_index = chunk.m_encoding_index;
      for (uint by = 0; by < 2; by++) {
        for (uint bx = 0; bx < 2; bx++) {
          if (g_tile_map[encoding_index][by][bx] == tile_index) {
            uint b = m_chunk_details[chunk_index].block_index[by][bx];
            uint64 selector = 0;
            for (uint sh = 0, p = 0; p < 16; p++, sh += 3) {
              uint8 s_best;
              for (uint32 error_best = UINT_MAX, t = 0; t < 8; t++) {
                uint8 s = alpha_order[t];
                int delta = m_blocks[b][p][component_index] - alpha_values[s];
                uint32 error = delta >= 0 ? delta : -delta;
                if (error < error_best) {
                  s_best = s;
                  error_best = error;
                }
              }
              selector |= (uint64)s_best << sh;
            }
            m_block_selectors[cAlpha0Chunks + alpha_index][b] = selector | encoding_weight[encoding_index] << 48;
            m_endpoint_indices[b].component[cAlpha0Chunks + alpha_index] = cluster_index;
          }
        }
      }
    }

    dxt_endpoint_refiner refiner;
    dxt_endpoint_refiner::params p;
    dxt_endpoint_refiner::results r;
    p.m_perceptual = m_params.m_perceptual;
    p.m_pSelectors = cluster.m_selectors.get_ptr();
    p.m_pPixels = cluster.m_pixels.get_ptr();
    p.m_num_pixels = cluster.m_pixels.size();
    p.m_dxt1_selectors = false;
    p.m_error_to_beat = cluster.m_error;
    p.m_block_index = cluster_index;
    cluster.m_refined.result = refiner.refine(p, r);
    cluster.m_refined.first_endpoint = r.m_low_color;
    cluster.m_refined.second_endpoint = r.m_high_color;
    cluster.m_refined.error = r.m_error;

    for (uint tile_iter = 0; tile_iter < cluster.m_tiles.size(); tile_iter++) {
      const uint chunk_index = cluster.m_tiles[tile_iter].first;
      const uint tile_index = cluster.m_tiles[tile_iter].second & 0xFFFFU;
      const uint alpha_index = cluster.m_tiles[tile_iter].second >> 16U;
      CRNLIB_ASSERT(chunk_index < m_num_chunks);
      CRNLIB_ASSERT(tile_index < cChunkMaxTiles);
      CRNLIB_ASSERT(alpha_index < m_num_alpha_blocks);

      compressed_chunk& chunk = m_compressed_chunks[cAlpha0Chunks + alpha_index][chunk_index];

      CRNLIB_ASSERT(chunk.m_endpoint_cluster_index[tile_index] == cluster_index);

      CRNLIB_ASSERT(tile_index < chunk.m_num_tiles);
      const compressed_tile& tile = chunk.m_tiles[tile_index];

      const chunk_tile_desc& layout = g_chunk_tile_layouts[tile.m_layout_index];
      layout;

      compressed_tile& quantized_tile = chunk.m_quantized_tiles[tile_index];

      quantized_tile.m_endpoint_cluster_index = cluster_index;
      quantized_tile.m_first_endpoint = results.m_first_endpoint;
      quantized_tile.m_second_endpoint = results.m_second_endpoint;
      quantized_tile.m_pixel_width = tile.m_pixel_width;
      quantized_tile.m_pixel_height = tile.m_pixel_height;
      quantized_tile.m_layout_index = tile.m_layout_index;
    }
  }
}

bool dxt_hc::determine_alpha_endpoint_codebook() {
  if (!m_num_alpha_blocks)
    return true;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Computing optimal alpha cluster endpoints");
#endif

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::determine_alpha_endpoint_codebook_task, i, NULL);

  m_pTask_pool->join();

  return !m_canceled;
}

void dxt_hc::create_selector_codebook_task(uint64 data, void* pData_ptr) {
  const uint thread_index = static_cast<uint>(data);
  const create_selector_codebook_state& state = *static_cast<create_selector_codebook_state*>(pData_ptr);

  for (uint comp_chunk_index = state.m_comp_index_start; comp_chunk_index <= state.m_comp_index_end; comp_chunk_index++) {
    const uint alpha_index = state.m_alpha_blocks ? (comp_chunk_index - cAlpha0Chunks) : 0;
    const uint alpha_pixel_comp = state.m_alpha_blocks ? m_params.m_alpha_component_indices[alpha_index] : 0;

    for (uint chunk_index = 0; chunk_index < m_num_chunks; chunk_index++) {
      if (m_canceled)
        return;

      if ((crn_get_current_thread_id() == m_main_thread_id) && ((chunk_index & 127) == 0)) {
        if (!update_progress(12 + comp_chunk_index, chunk_index, m_num_chunks))
          return;
      }

      if (m_pTask_pool->get_num_threads()) {
        if ((chunk_index % (m_pTask_pool->get_num_threads() + 1)) != thread_index)
          continue;
      }

      compressed_chunk& chunk = m_compressed_chunks[comp_chunk_index][chunk_index];

      for (uint tile_index = 0; tile_index < chunk.m_num_tiles; tile_index++) {
        compressed_tile& quantized_tile = chunk.m_quantized_tiles[tile_index];

        const chunk_tile_desc& layout = g_chunk_tile_layouts[quantized_tile.m_layout_index];

        const uint tile_blocks_x = layout.m_width >> 2;
        const uint tile_blocks_y = layout.m_height >> 2;

        const uint tile_block_ofs_x = layout.m_x_ofs >> 2;
        const uint tile_block_ofs_y = layout.m_y_ofs >> 2;

        if (state.m_alpha_blocks) {
          uint block_values[cDXT5SelectorValues];
          dxt5_block::get_block_values(block_values, quantized_tile.m_first_endpoint, quantized_tile.m_second_endpoint);

          for (uint by = 0; by < tile_blocks_y; by++) {
            for (uint bx = 0; bx < tile_blocks_x; bx++) {
#if 0
                        uint best_index = selector_vq.find_best_codebook_entry_fs(training_vecs[comp_chunk_index][(tile_block_ofs_x+bx)+(tile_block_ofs_y+by)*2][chunk_index]);
#else
              const dxt_pixel_block& block = m_pChunks[chunk_index].m_blocks[tile_block_ofs_y + by][tile_block_ofs_x + bx];

              uint best_error = UINT_MAX;
              uint best_index = 0;

              for (uint i = 0; i < state.m_selectors_cb.size(); i++) {
                const selectors& s = state.m_selectors_cb[i];

                uint total_error = 0;

                for (uint y = 0; y < cBlockPixelHeight; y++) {
                  for (uint x = 0; x < cBlockPixelWidth; x++) {
                    int a = block.m_pixels[y][x][alpha_pixel_comp];
                    int b = block_values[s.m_selectors[y][x]];
                    int error = a - b;
                    error *= error;

                    total_error += error;
                    if (total_error > best_error)
                      goto early_out;
                  }  // x
                }    //y

              early_out:
                if (total_error < best_error) {
                  best_error = total_error;
                  best_index = i;

                  if (best_error == 0)
                    break;
                }
              }  // i
#endif

              CRNLIB_ASSERT((tile_block_ofs_x + bx) < 2);
              CRNLIB_ASSERT((tile_block_ofs_y + by) < 2);

              chunk.m_selector_cluster_index[tile_block_ofs_y + by][tile_block_ofs_x + bx] = static_cast<uint16>(best_index);

              {
                scoped_spinlock lock(state.m_chunk_blocks_using_selectors_lock);
                state.m_chunk_blocks_using_selectors[best_index].push_back(block_id(chunk_index, alpha_index, tile_index, tile_block_ofs_x + bx, tile_block_ofs_y + by));
              }
              //   std::make_pair(chunk_index, (tile_index << 16) | ((tile_block_ofs_y + by) << 8) | (tile_block_ofs_x + bx) ) );

            }  // bx
          }    // by

        } else {
          color_quad_u8 block_colors[cDXT1SelectorValues];
          dxt1_block::get_block_colors4(block_colors, static_cast<uint16>(quantized_tile.m_first_endpoint), static_cast<uint16>(quantized_tile.m_second_endpoint));

          const bool block_with_alpha = quantized_tile.m_first_endpoint == quantized_tile.m_second_endpoint;

          for (uint by = 0; by < tile_blocks_y; by++) {
            for (uint bx = 0; bx < tile_blocks_x; bx++) {
              const dxt_pixel_block& block = m_pChunks[chunk_index].m_blocks[tile_block_ofs_y + by][tile_block_ofs_x + bx];

              uint best_error = UINT_MAX;
              uint best_index = 0;

              for (uint i = 0; i < state.m_selectors_cb.size(); i++) {
                const selectors& s = state.m_selectors_cb[i];

                uint total_error = 0;

                for (uint y = 0; y < cBlockPixelHeight; y++) {
                  for (uint x = 0; x < cBlockPixelWidth; x++) {
                    const color_quad_u8& a = block.m_pixels[y][x];

                    uint selector_index = s.m_selectors[y][x];
                    if ((block_with_alpha) && (selector_index == 3))
                      total_error += 999999;

                    const color_quad_u8& b = block_colors[selector_index];

                    uint error = color::color_distance(m_params.m_perceptual, a, b, false);

                    total_error += error;
                    if (total_error > best_error)
                      goto early_out2;
                  }  // x
                }    //y

              early_out2:
                if (total_error < best_error) {
                  best_error = total_error;
                  best_index = i;

                  if (best_error == 0)
                    break;
                }
              }  // i

              CRNLIB_ASSERT((tile_block_ofs_x + bx) < 2);
              CRNLIB_ASSERT((tile_block_ofs_y + by) < 2);

              chunk.m_selector_cluster_index[tile_block_ofs_y + by][tile_block_ofs_x + bx] = static_cast<uint16>(best_index);

              {
                scoped_spinlock lock(state.m_chunk_blocks_using_selectors_lock);
                state.m_chunk_blocks_using_selectors[best_index].push_back(block_id(chunk_index, 0, tile_index, tile_block_ofs_x + bx, tile_block_ofs_y + by));
              }
              //   std::make_pair(chunk_index, (tile_index << 16) | ((tile_block_ofs_y + by) << 8) | (tile_block_ofs_x + bx) ) );

            }  // bx
          }    // by

        }  // if alpha_blocks

      }  // tile_index

    }  // chunk_index

  }  // comp_chunk_index
}

bool dxt_hc::create_selector_codebook(bool alpha_blocks) {
  vec16F_tree_vq selector_vq;
  vec16F v;
  uint c_start = alpha_blocks ? cAlpha0Chunks : cColorChunks;
  uint c_end = alpha_blocks ? cAlpha0Chunks + m_num_alpha_blocks - 1 : cColorChunks;
  float scale = alpha_blocks ? 0.125f : 0.25f;

  for (uint c = c_start; c <= c_end; c++) {
    for (uint b = 0; b < m_blocks.size(); b++) {
      uint64 selector = m_block_selectors[c][b];
      for (uint8 p = 0; p < 16; p++, selector >>= 3)
        v[p] = ((selector & 7) + 0.5f) * scale;
      selector_vq.add_training_vec(v, selector);
    }
  }

  selector_vq.generate_codebook(alpha_blocks ? m_params.m_alpha_selector_codebook_size : m_params.m_color_selector_codebook_size);
  selectors_vec& selectors_cb = alpha_blocks ? m_alpha_selectors : m_color_selectors;
  selectors_cb.resize(selector_vq.get_codebook_size());

  for (uint i = 0; i < selector_vq.get_codebook_size(); i++) {
    const vec16F& v = selector_vq.get_codebook_entry(i);
    for (uint j = 0; j < 16; j++)
      selectors_cb[i].m_selectors[j >> 2][j & 3] = alpha_blocks ? g_dxt5_from_linear[(int)(v[j] * 8.0f)] : g_dxt1_from_linear[(int)(v[j] * 4.0f)];
  }

  chunk_blocks_using_selectors_vec& chunk_blocks_using_selectors = alpha_blocks ? m_chunk_blocks_using_alpha_selectors : m_chunk_blocks_using_color_selectors;

  chunk_blocks_using_selectors.clear();
  chunk_blocks_using_selectors.resize(selectors_cb.size());

  create_selector_codebook_state state(*this, alpha_blocks, c_start, c_end, selector_vq, chunk_blocks_using_selectors, selectors_cb);

  for (uint i = 0; i <= m_pTask_pool->get_num_threads(); i++)
    m_pTask_pool->queue_object_task(this, &dxt_hc::create_selector_codebook_task, i, &state);

  m_pTask_pool->join();

  return !m_canceled;
}

bool dxt_hc::refine_quantized_color_selectors() {
  if (!m_has_color_blocks)
    return true;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Refining quantized color selectors");
#endif

  uint total_refined_selectors = 0;
  uint total_refined_pixels = 0;
  uint total_selectors = 0;

  for (uint selector_index = 0; selector_index < m_color_selectors.size(); selector_index++) {
    if ((selector_index & 255) == 0) {
      if (!update_progress(15, selector_index, m_color_selectors.size()))
        return false;
    }

    if (m_chunk_blocks_using_color_selectors[selector_index].empty())
      continue;

    selectors& sel = m_color_selectors[selector_index];

    for (uint y = 0; y < cBlockPixelHeight; y++) {
      for (uint x = 0; x < cBlockPixelWidth; x++) {
        uint best_s = 0;
        uint best_error = UINT_MAX;

        for (uint s = 0; s < cDXT1SelectorValues; s++) {
          uint total_error = 0;

          for (uint block_iter = 0; block_iter < m_chunk_blocks_using_color_selectors[selector_index].size(); block_iter++) {
            const block_id& id = m_chunk_blocks_using_color_selectors[selector_index][block_iter];
            const uint chunk_index = id.m_chunk_index;
            const uint tile_index = id.m_tile_index;
            const uint chunk_block_x = id.m_block_x;
            const uint chunk_block_y = id.m_block_y;

            CRNLIB_ASSERT((chunk_block_x < cChunkBlockWidth) && (chunk_block_y < cChunkBlockHeight));

            const compressed_chunk& chunk = m_compressed_chunks[cColorChunks][chunk_index];
            CRNLIB_ASSERT(tile_index < chunk.m_num_tiles);

            CRNLIB_ASSERT(chunk.m_selector_cluster_index[chunk_block_y][chunk_block_x] == selector_index);

            const compressed_tile& tile = chunk.m_quantized_tiles[tile_index];

            //const chunk_tile_desc& tile_desc = g_chunk_tile_layouts[tile.m_layout_index];

            color_quad_u8 block_colors[cDXT1SelectorValues];
            CRNLIB_ASSERT(tile.m_first_endpoint >= tile.m_second_endpoint);
            dxt1_block::get_block_colors4(block_colors, static_cast<uint16>(tile.m_first_endpoint), static_cast<uint16>(tile.m_second_endpoint));

            if ((tile.m_first_endpoint == tile.m_second_endpoint) && (s == 3))
              total_error += 999999;

            const color_quad_u8& orig_pixel = m_pChunks[chunk_index](chunk_block_x * cBlockPixelWidth + x, chunk_block_y * cBlockPixelHeight + y);
            const color_quad_u8& quantized_pixel = block_colors[s];

            const uint error = color::color_distance(m_params.m_perceptual, orig_pixel, quantized_pixel, false);
            total_error += error;

          }  // block_iter

          if (total_error < best_error) {
            best_error = total_error;
            best_s = s;
          }

        }  // s

        if (sel.m_selectors[y][x] != best_s) {
          total_refined_selectors++;
          total_refined_pixels += m_chunk_blocks_using_color_selectors[selector_index].size();
          sel.m_selectors[y][x] = static_cast<uint8>(best_s);
        }

        total_selectors++;

      }  //x

    }  //y

  }  // selector_index

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Total refined pixels: %u, selectors: %u out of %u", total_refined_pixels, total_refined_selectors, total_selectors);
#endif

  uint selector_count = 0;
  hash_map<uint32, uint> packed_selectors;
  for (uint i = 0; i < m_color_selectors.size(); i++) {
    if (m_chunk_blocks_using_color_selectors[i].size()) {
      uint32 packed_selector = 0;
      for (uint s = 0; s < 16; s++)
        packed_selector |= m_color_selectors[i].get_by_index(s) << (s << 1);
      hash_map<uint32, uint>::insert_result insert_result = packed_selectors.insert(packed_selector, selector_count);
      if (insert_result.second) {
        if (selector_count != i) {
          m_color_selectors[selector_count] = m_color_selectors[i];
          m_chunk_blocks_using_color_selectors[selector_count].swap(m_chunk_blocks_using_color_selectors[i]);
          for (uint b = 0; b < m_chunk_blocks_using_color_selectors[selector_count].size(); b++) {
            const block_id& id = m_chunk_blocks_using_color_selectors[selector_count][b];
            m_compressed_chunks[cColorChunks][id.m_chunk_index].m_selector_cluster_index[id.m_block_y][id.m_block_x] = selector_count;
          }
        }
        selector_count++;
      } else {
        for (uint b = 0; b < m_chunk_blocks_using_color_selectors[i].size(); b++) {
          const block_id& id = m_chunk_blocks_using_color_selectors[i][b];
          m_compressed_chunks[cColorChunks][id.m_chunk_index].m_selector_cluster_index[id.m_block_y][id.m_block_x] = insert_result.first->second;
        }
        m_chunk_blocks_using_color_selectors[insert_result.first->second].append(m_chunk_blocks_using_color_selectors[i]);
        m_chunk_blocks_using_color_selectors[i].clear();
      }
    }
  }
  m_color_selectors.resize(selector_count);
  m_chunk_blocks_using_color_selectors.resize(selector_count);
  
  return true;
}

bool dxt_hc::refine_quantized_alpha_selectors() {
  if (!m_num_alpha_blocks)
    return true;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Refining quantized alpha selectors");
#endif

  uint total_refined_selectors = 0;
  uint total_refined_pixels = 0;
  uint total_selectors = 0;

  for (uint selector_index = 0; selector_index < m_alpha_selectors.size(); selector_index++) {
    if ((selector_index & 255) == 0) {
      if (!update_progress(16, selector_index, m_alpha_selectors.size()))
        return false;
    }

    if (m_chunk_blocks_using_alpha_selectors[selector_index].empty())
      continue;

    selectors& sel = m_alpha_selectors[selector_index];

    for (uint y = 0; y < cBlockPixelHeight; y++) {
      for (uint x = 0; x < cBlockPixelWidth; x++) {
        uint best_s = 0;
        uint best_error = UINT_MAX;

        for (uint s = 0; s < cDXT5SelectorValues; s++) {
          uint total_error = 0;

          for (uint block_iter = 0; block_iter < m_chunk_blocks_using_alpha_selectors[selector_index].size(); block_iter++) {
            const block_id& id = m_chunk_blocks_using_alpha_selectors[selector_index][block_iter];
            const uint chunk_index = id.m_chunk_index;
            const uint tile_index = id.m_tile_index;
            const uint chunk_block_x = id.m_block_x;
            const uint chunk_block_y = id.m_block_y;
            const uint alpha_index = id.m_alpha_index;
            CRNLIB_ASSERT(alpha_index < m_num_alpha_blocks);

            CRNLIB_ASSERT((chunk_block_x < cChunkBlockWidth) && (chunk_block_y < cChunkBlockHeight));

            const compressed_chunk& chunk = m_compressed_chunks[alpha_index + cAlpha0Chunks][chunk_index];
            CRNLIB_ASSERT(tile_index < chunk.m_num_tiles);

            CRNLIB_ASSERT(chunk.m_selector_cluster_index[chunk_block_y][chunk_block_x] == selector_index);

            const compressed_tile& tile = chunk.m_quantized_tiles[tile_index];

            //const chunk_tile_desc& tile_desc = g_chunk_tile_layouts[tile.m_layout_index];

            uint block_values[cDXT5SelectorValues];
            CRNLIB_ASSERT(tile.m_first_endpoint >= tile.m_second_endpoint);
            dxt5_block::get_block_values(block_values, tile.m_first_endpoint, tile.m_second_endpoint);

            int orig_value = m_pChunks[chunk_index](chunk_block_x * cBlockPixelWidth + x, chunk_block_y * cBlockPixelHeight + y)[m_params.m_alpha_component_indices[alpha_index]];
            int quantized_value = block_values[s];

            int error = (orig_value - quantized_value);
            error *= error;

            total_error += error;

          }  // block_iter

          if (total_error < best_error) {
            best_error = total_error;
            best_s = s;
          }

        }  // s

        if (sel.m_selectors[y][x] != best_s) {
          total_refined_selectors++;
          total_refined_pixels += m_chunk_blocks_using_alpha_selectors[selector_index].size();
          sel.m_selectors[y][x] = static_cast<uint8>(best_s);
        }

        total_selectors++;

      }  //x

    }  //y

  }  // selector_index

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    console::info("Total refined pixels: %u, selectors: %u out of %u", total_refined_pixels, total_refined_selectors, total_selectors);
#endif

  uint selector_count = 0;
  hash_map<uint64, uint> packed_selectors;
  for (uint i = 0; i < m_alpha_selectors.size(); i++) {
    if (m_chunk_blocks_using_alpha_selectors[i].size()) {
      uint64 packed_selector = 0;
      for (uint s = 0; s < 16; s++)
        packed_selector |= (uint64)m_alpha_selectors[i].get_by_index(s) << (s << 2);
      hash_map<uint64, uint>::insert_result insert_result = packed_selectors.insert(packed_selector, selector_count);
      if (insert_result.second) {
        if (selector_count != i) {
          m_alpha_selectors[selector_count] = m_alpha_selectors[i];
          m_chunk_blocks_using_alpha_selectors[selector_count].swap(m_chunk_blocks_using_alpha_selectors[i]);
          for (uint b = 0; b < m_chunk_blocks_using_alpha_selectors[selector_count].size(); b++) {
            const block_id& id = m_chunk_blocks_using_alpha_selectors[selector_count][b];
            m_compressed_chunks[cAlpha0Chunks + id.m_alpha_index][id.m_chunk_index].m_selector_cluster_index[id.m_block_y][id.m_block_x] = selector_count;
          }
        }
        selector_count++;
      } else {
        for (uint b = 0; b < m_chunk_blocks_using_alpha_selectors[i].size(); b++) {
          const block_id& id = m_chunk_blocks_using_alpha_selectors[i][b];
          m_compressed_chunks[cAlpha0Chunks + id.m_alpha_index][id.m_chunk_index].m_selector_cluster_index[id.m_block_y][id.m_block_x] = insert_result.first->second;
        }
        m_chunk_blocks_using_alpha_selectors[insert_result.first->second].append(m_chunk_blocks_using_alpha_selectors[i]);
        m_chunk_blocks_using_alpha_selectors[i].clear();
      }
    }
  }
  m_alpha_selectors.resize(selector_count);
  m_chunk_blocks_using_alpha_selectors.resize(selector_count);
  
  return true;
}

bool dxt_hc::refine_quantized_color_endpoints() {
  if (!m_has_color_blocks)
    return true;

  for (uint cluster_index = 0; cluster_index < m_color_clusters.size(); cluster_index++) {
    tile_cluster& cluster = m_color_clusters[cluster_index];
    if (cluster.m_refined.result) {
      cluster.m_error = cluster.m_refined.error;
      cluster.m_first_endpoint = cluster.m_refined.first_endpoint;
      cluster.m_second_endpoint = cluster.m_refined.second_endpoint;
      for (uint tile_iter = 0; tile_iter < cluster.m_tiles.size(); tile_iter++) {
        const uint chunk_index = cluster.m_tiles[tile_iter].first;
        const uint tile_index = cluster.m_tiles[tile_iter].second;
        compressed_chunk& chunk = m_compressed_chunks[cColorChunks][chunk_index];
        compressed_tile& tile = chunk.m_quantized_tiles[tile_index];
        tile.m_first_endpoint = cluster.m_first_endpoint;
        tile.m_second_endpoint = cluster.m_second_endpoint;
      }
    }
  }  

  uint cluster_count = 0;
  hash_map<uint32, uint> packed_clusters;
  for (uint i = 0; i < m_color_clusters.size(); i++) {
    tile_cluster& cluster = m_color_clusters[i];
    if (cluster.m_tiles.size()) {
      uint32 packed_cluster = cluster.m_first_endpoint | cluster.m_second_endpoint << 16;
      hash_map<uint32, uint>::insert_result insert_result = packed_clusters.insert(packed_cluster, cluster_count);
      if (insert_result.second) {
        if (cluster_count != i) {
          tile_cluster& destination_cluster = m_color_clusters[cluster_count];
          destination_cluster.m_error = cluster.m_error;
          destination_cluster.m_first_endpoint = cluster.m_first_endpoint;
          destination_cluster.m_second_endpoint = cluster.m_second_endpoint;
          destination_cluster.m_tiles.swap(cluster.m_tiles);
          for (uint t = 0; t < destination_cluster.m_tiles.size(); t++) {
            const uint chunk_index = destination_cluster.m_tiles[t].first;
            const uint tile_index = destination_cluster.m_tiles[t].second;
            compressed_tile& tile = m_compressed_chunks[cColorChunks][chunk_index].m_quantized_tiles[tile_index];
            tile.m_first_endpoint = destination_cluster.m_first_endpoint;
            tile.m_second_endpoint = destination_cluster.m_second_endpoint;
            tile.m_endpoint_cluster_index = cluster_count;
          }
        }
        cluster_count++;
      } else {
        tile_cluster& destination_cluster = m_color_clusters[insert_result.first->second];
        for (uint t = 0; t < cluster.m_tiles.size(); t++) {
          const uint chunk_index = cluster.m_tiles[t].first;
          const uint tile_index = cluster.m_tiles[t].second;
          compressed_tile& tile = m_compressed_chunks[cColorChunks][chunk_index].m_quantized_tiles[tile_index];
          tile.m_endpoint_cluster_index = insert_result.first->second;
        }
        destination_cluster.m_error += cluster.m_error;
        destination_cluster.m_tiles.append(cluster.m_tiles);
        cluster.m_tiles.clear();
      }
    }
  }
  m_color_clusters.resize(cluster_count);

  return true;
}

bool dxt_hc::refine_quantized_alpha_endpoints() {
  if (!m_num_alpha_blocks)
    return true;

  for (uint cluster_index = 0; cluster_index < m_alpha_clusters.size(); cluster_index++) {
    tile_cluster& cluster = m_alpha_clusters[cluster_index];
    if (cluster.m_refined.result) {
      cluster.m_error = cluster.m_refined.error;
      cluster.m_first_endpoint = cluster.m_refined.first_endpoint;
      cluster.m_second_endpoint = cluster.m_refined.second_endpoint;
      for (uint tile_iter = 0; tile_iter < cluster.m_tiles.size(); tile_iter++) {
        const uint chunk_index = cluster.m_tiles[tile_iter].first;
        const uint tile_index = cluster.m_tiles[tile_iter].second & 0xFFFFU;
        const uint alpha_index = cluster.m_tiles[tile_iter].second >> 16U;
        compressed_tile& tile = m_compressed_chunks[cAlpha0Chunks + alpha_index][chunk_index].m_quantized_tiles[tile_index];
        tile.m_first_endpoint = cluster.m_first_endpoint;
        tile.m_second_endpoint = cluster.m_second_endpoint;
      }
    }
  }

  uint cluster_count = 0;
  hash_map<uint32, uint> packed_clusters;
  for (uint i = 0; i < m_alpha_clusters.size(); i++) {
    tile_cluster& cluster = m_alpha_clusters[i];
    if (cluster.m_tiles.size()) {
      uint32 packed_cluster = cluster.m_first_endpoint | cluster.m_second_endpoint << 16;
      hash_map<uint32, uint>::insert_result insert_result = packed_clusters.insert(packed_cluster, cluster_count);
      if (insert_result.second) {
        if (cluster_count != i) {
          tile_cluster& destination_cluster = m_alpha_clusters[cluster_count];
          destination_cluster.m_error = cluster.m_error;
          destination_cluster.m_first_endpoint = cluster.m_first_endpoint;
          destination_cluster.m_second_endpoint = cluster.m_second_endpoint;
          destination_cluster.m_tiles.swap(cluster.m_tiles);
          for (uint t = 0; t < destination_cluster.m_tiles.size(); t++) {
            const uint chunk_index = destination_cluster.m_tiles[t].first;
            const uint tile_index = destination_cluster.m_tiles[t].second & 0xFFFFU;
            const uint alpha_index = destination_cluster.m_tiles[t].second >> 16U;
            compressed_tile& tile = m_compressed_chunks[cAlpha0Chunks + alpha_index][chunk_index].m_quantized_tiles[tile_index];
            tile.m_first_endpoint = destination_cluster.m_first_endpoint;
            tile.m_second_endpoint = destination_cluster.m_second_endpoint;
            tile.m_endpoint_cluster_index = cluster_count;
          }
        }
        cluster_count++;
      } else {
        tile_cluster& destination_cluster = m_alpha_clusters[insert_result.first->second];
        for (uint t = 0; t < cluster.m_tiles.size(); t++) {
          const uint chunk_index = cluster.m_tiles[t].first;
          const uint tile_index = cluster.m_tiles[t].second & 0xFFFFU;
          const uint alpha_index = cluster.m_tiles[t].second >> 16U;
          compressed_tile& tile = m_compressed_chunks[cAlpha0Chunks + alpha_index][chunk_index].m_quantized_tiles[tile_index];
          tile.m_endpoint_cluster_index = insert_result.first->second;
        }
        destination_cluster.m_error += cluster.m_error;
        destination_cluster.m_tiles.append(cluster.m_tiles);
        cluster.m_tiles.clear();
      }
    }
  }
  m_alpha_clusters.resize(cluster_count);

  return true;
}

bool dxt_hc::create_block_encodings(const params& p) {
  crnlib::vector<endpoint_indices_details>& endpoint_indices = *p.m_endpoint_indices;
  crnlib::vector<selector_indices_details>& selector_indices = *p.m_selector_indices;

  endpoint_indices.resize(m_num_chunks << 2);
  selector_indices.resize(m_num_chunks << 2);
  bool hasBlocks[cNumCompressedChunkVecs] = {m_has_color_blocks, m_num_alpha_blocks > 0, m_num_alpha_blocks > 1};

  for (uint level = 0; level < p.m_num_levels; level++) {
    uint first_chunk = p.m_levels[level].m_first_chunk;
    uint end_chunk = p.m_levels[level].m_first_chunk + p.m_levels[level].m_num_chunks;
    uint chunk_width = p.m_levels[level].m_chunk_width;
    uint block_width = chunk_width << 1;
    for (uint b = first_chunk << 2, cy = 0, chunk_base = first_chunk; chunk_base < end_chunk; chunk_base += chunk_width, cy++) {
      for (uint by = 0; by < 2; by++) {
        for (uint cx = 0; cx < chunk_width; cx++) {
          for (uint bx = 0; bx < 2; bx++, b++) {
            bool top_match = cy || by;
            bool left_match = top_match || cx || bx;
            for (uint c = 0; c < cNumCompressedChunkVecs; c++) {
              if (hasBlocks[c]) {
                const compressed_chunk& chunk = m_compressed_chunks[c][chunk_base + cx];
                uint16 endpoint_index = chunk.m_quantized_tiles[g_tile_map[chunk.m_encoding_index][by][bx]].m_endpoint_cluster_index;
                left_match = left_match && endpoint_index == endpoint_indices[b - 1].component[c];
                top_match = top_match && endpoint_index == endpoint_indices[b - block_width].component[c];
                endpoint_indices[b].component[c] = endpoint_index;
                selector_indices[b].component[c] = chunk.m_selector_cluster_index[by][bx];
              }
            }
            endpoint_indices[b].reference = left_match ? 1 : top_match ? 2 : 0;
          }
        }
      }
    }
  }

  if (m_has_color_blocks) {
    m_color_endpoints.resize(m_color_clusters.size());
    for (uint i = 0; i < m_color_clusters.size(); i++)
      m_color_endpoints[i] = dxt1_block::pack_endpoints(m_color_clusters[i].m_first_endpoint, m_color_clusters[i].m_second_endpoint);
  }

  if (m_num_alpha_blocks) {
    m_alpha_endpoints.resize(m_alpha_clusters.size());
    for (uint i = 0; i < m_alpha_clusters.size(); i++)
      m_alpha_endpoints[i] = dxt5_block::pack_endpoints(m_alpha_clusters[i].m_first_endpoint, m_alpha_clusters[i].m_second_endpoint);
  }

  return true;
}

bool dxt_hc::update_progress(uint phase_index, uint subphase_index, uint subphase_total) {
  CRNLIB_ASSERT(crn_get_current_thread_id() == m_main_thread_id);

  if (!m_params.m_pProgress_func)
    return true;

#if CRNLIB_ENABLE_DEBUG_MESSAGES
  if (m_params.m_debugging)
    return true;
#endif

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
