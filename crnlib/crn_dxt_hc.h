// File: crn_dxt_hc.h
// See Copyright Notice and license at the end of inc/crnlib.h
#pragma once
#include "crn_dxt1.h"
#include "crn_dxt5a.h"
#include "crn_dxt_endpoint_refiner.h"
#include "crn_image.h"
#include "crn_dxt.h"
#include "crn_image.h"
#include "crn_dxt_hc_common.h"
#include "crn_tree_clusterizer.h"
#include "crn_threading.h"

#define CRN_NO_FUNCTION_DEFINITIONS
#include "../inc/crnlib.h"

namespace crnlib {
const uint cTotalCompressionPhases = 25;

class dxt_hc {
 public:
  dxt_hc();
  ~dxt_hc();

  struct endpoint_indices_details {
    union {
      struct {
        uint16 color;
        uint16 alpha0;
        uint16 alpha1;
      };
      uint16 component[3];
    };
    uint8 reference;
  };

  struct selector_indices_details {
    union {
      struct {
        uint16 color;
        uint16 alpha0;
        uint16 alpha1;
      };
      uint16 component[3];
    };
  };

  struct tile_details {
    crnlib::vector<color_quad_u8> pixels;
    float weight;
    vec<6, float> color_endpoint;
    vec<2, float> alpha_endpoints[2];
    uint16 cluster_indices[3];
  };
  crnlib::vector<tile_details> m_tiles;
  uint m_num_tiles;
  float m_color_derating[cCRNMaxLevels][8];
  float m_alpha_derating[8];

  color_quad_u8 (*m_blocks)[16];
  uint m_num_blocks;
  crnlib::vector<float> m_block_weights;
  crnlib::vector<uint8> m_block_encodings;
  crnlib::vector<uint64> m_block_selectors[3];
  crnlib::vector<uint32> m_color_selectors;
  crnlib::vector<uint64> m_alpha_selectors;
  crnlib::vector<bool> m_color_selectors_used;
  crnlib::vector<bool> m_alpha_selectors_used;
  crnlib::vector<uint> m_tile_indices;
  crnlib::vector<endpoint_indices_details> m_endpoint_indices;
  crnlib::vector<selector_indices_details> m_selector_indices;

  struct params {
    params()
        : m_blocks(0),
          m_num_blocks(0),
          m_num_levels(0),
          m_num_faces(0),
          m_format(cDXT1),
          m_perceptual(true),
          m_hierarchical(true),
          m_color_endpoint_codebook_size(3072),
          m_color_selector_codebook_size(3072),
          m_alpha_endpoint_codebook_size(3072),
          m_alpha_selector_codebook_size(3072),
          m_adaptive_tile_color_psnr_derating(2.0f),
          m_adaptive_tile_alpha_psnr_derating(2.0f),
          m_adaptive_tile_color_alpha_weighting_ratio(3.0f),
          m_debugging(false),
          m_pProgress_func(0),
          m_pProgress_func_data(0),
          m_endpoint_indices(0),
          m_selector_indices(0) {
      m_alpha_component_indices[0] = 3;
      m_alpha_component_indices[1] = 0;
      for (uint i = 0; i < cCRNMaxLevels; i++) {
        m_levels[i].m_first_block = 0;
        m_levels[i].m_num_blocks = 0;
        m_levels[i].m_block_width = 0;
      }
    }

    color_quad_u8 (*m_blocks)[16];
    uint m_num_blocks;
    uint m_num_levels;
    uint m_num_faces;

    struct {
      uint m_first_block;
      uint m_num_blocks;
      uint m_block_width;
      float m_weight;
    } m_levels[cCRNMaxLevels];

    dxt_format m_format;
    bool m_perceptual;
    bool m_hierarchical;

    uint m_color_endpoint_codebook_size;
    uint m_color_selector_codebook_size;
    uint m_alpha_endpoint_codebook_size;
    uint m_alpha_selector_codebook_size;

    float m_adaptive_tile_color_psnr_derating;
    float m_adaptive_tile_alpha_psnr_derating;
    float m_adaptive_tile_color_alpha_weighting_ratio;
    uint m_alpha_component_indices[2];

    bool m_debugging;
    crn_progress_callback_func m_pProgress_func;
    void* m_pProgress_func_data;

    crnlib::vector<endpoint_indices_details> *m_endpoint_indices;
    crnlib::vector<selector_indices_details> *m_selector_indices;
  };

  struct selectors {
    selectors() { utils::zero_object(*this); }

    uint8 m_selectors[cBlockPixelHeight][cBlockPixelWidth];

    uint8 get_by_index(uint i) const {
      CRNLIB_ASSERT(i < (cBlockPixelWidth * cBlockPixelHeight));
      const uint8* p = (const uint8*)m_selectors;
      return *(p + i);
    }
    void set_by_index(uint i, uint v) {
      CRNLIB_ASSERT(i < (cBlockPixelWidth * cBlockPixelHeight));
      uint8* p = (uint8*)m_selectors;
      *(p + i) = static_cast<uint8>(v);
    }
  };
  typedef crnlib::vector<selectors> selectors_vec;

  void clear();
  bool compress(const params& p, task_pool& task_pool);

  // Color endpoints
  inline uint get_color_endpoint_codebook_size() const { return m_color_endpoints.size(); }
  inline uint get_color_endpoint(uint codebook_index) const { return m_color_endpoints[codebook_index]; }
  const crnlib::vector<uint>& get_color_endpoint_vec() const { return m_color_endpoints; }

  // Color selectors
  uint get_color_selector_codebook_size() const { return m_color_selectors_vec.size(); }
  const selectors& get_color_selectors(uint codebook_index) const { return m_color_selectors_vec[codebook_index]; }
  const crnlib::vector<selectors>& get_color_selectors_vec() const { return m_color_selectors_vec; }

  // Alpha endpoints
  inline uint get_alpha_endpoint_codebook_size() const { return m_alpha_endpoints.size(); }
  inline uint get_alpha_endpoint(uint codebook_index) const { return m_alpha_endpoints[codebook_index]; }
  const crnlib::vector<uint>& get_alpha_endpoint_vec() const { return m_alpha_endpoints; }

  // Alpha selectors
  uint get_alpha_selector_codebook_size() const { return m_alpha_selectors_vec.size(); }
  const selectors& get_alpha_selectors(uint codebook_index) const { return m_alpha_selectors_vec[codebook_index]; }
  const crnlib::vector<selectors>& get_alpha_selectors_vec() const { return m_alpha_selectors_vec; }

 private:
  params m_params;

  uint m_num_alpha_blocks;
  bool m_has_color_blocks;
  bool m_has_alpha0_blocks;
  bool m_has_alpha1_blocks;

  enum {
    cColorBlocks = 0,
    cAlpha0Blocks = 1,
    cAlpha1Blocks = 2,
    cNumCompressedComponents = 3
  };

  struct endpoint_cluster {
    endpoint_cluster() : m_first_endpoint(0), m_second_endpoint(0) {}
    crnlib::vector<uint> m_blocks[3];
    crnlib::vector<color_quad_u8> m_pixels;
    uint m_first_endpoint;
    uint m_second_endpoint;
    color_quad_u8 m_color_values[4];
    uint m_alpha_values[8];
    bool m_refined_result;
    uint m_refined_first_endpoint;
    uint m_refined_second_endpoint;
    uint m_refined_alpha_values[8];
  };
  crnlib::vector<endpoint_cluster> m_color_clusters;
  crnlib::vector<endpoint_cluster> m_alpha_clusters;

  selectors_vec m_alpha_selectors_vec;
  selectors_vec m_color_selectors_vec;
  crnlib::vector<uint> m_color_endpoints;
  crnlib::vector<uint> m_alpha_endpoints;

  crn_thread_id_t m_main_thread_id;
  bool m_canceled;
  task_pool* m_pTask_pool;

  int m_prev_phase_index;
  int m_prev_percentage_complete;

  typedef vec<6, float> vec6F;
  typedef vec<16, float> vec16F;
  typedef tree_clusterizer<vec2F> vec2F_tree_vq;
  typedef tree_clusterizer<vec6F> vec6F_tree_vq;
  typedef tree_clusterizer<vec16F> vec16F_tree_vq;

  void determine_tiles_task(uint64 data, void* pData_ptr);

  void determine_color_endpoint_codebook_task(uint64 data, void* pData_ptr);
  void determine_color_endpoint_clusters_task(uint64 data, void* pData_ptr);
  void determine_color_endpoints();

  void determine_alpha_endpoint_codebook_task(uint64 data, void* pData_ptr);
  void determine_alpha_endpoint_clusters_task(uint64 data, void* pData_ptr);
  void determine_alpha_endpoints();

  void create_color_selector_codebook_task(uint64 data, void* pData_ptr);
  void create_color_selector_codebook();

  void create_alpha_selector_codebook_task(uint64 data, void* pData_ptr);
  void create_alpha_selector_codebook();

  bool update_progress(uint phase_index, uint subphase_index, uint subphase_total);
};

CRNLIB_DEFINE_BITWISE_COPYABLE(dxt_hc::selectors);

}  // namespace crnlib
