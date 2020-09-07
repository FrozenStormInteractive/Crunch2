configure_file(${CMAKE_CURRENT_SOURCE_DIR}/crn_version.h.in ${CMAKE_CURRENT_BINARY_DIR}/crn_version.h @ONLY)

set(CRNLIB_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_arealist.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_arealist.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_assert.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_assert.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_buffer_stream.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_cfile_stream.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_checksum.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_checksum.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_clusterizer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_color.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_colorized_console.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_colorized_console.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_command_line_params.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_command_line_params.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_comp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_comp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_console.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_console.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_core.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_core.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_data_stream.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_data_stream.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_data_stream_serializer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dds_comp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dds_comp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_decomp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt1.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt5a.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt5a.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt_endpoint_refiner.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt_endpoint_refiner.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt_fast.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt_fast.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt_hc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt_hc.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt_hc_common.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt_hc_common.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt_image.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dxt_image.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dynamic_stream.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dynamic_string.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_dynamic_string.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_etc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_etc.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_file_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_file_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_find_files.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_find_files.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_hash.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_hash.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_hash_map.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_hash_map.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_helpers.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_huffman_codes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_huffman_codes.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_image.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_image_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_image_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_intersect.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_ktx_texture.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_ktx_texture.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_lzma_codec.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_lzma_codec.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_math.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_math.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_matrix.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_mem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_mem.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_mipmapped_texture.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_mipmapped_texture.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_packed_uint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_pixel_format.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_pixel_format.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_platform.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_platform.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_prefix_coding.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_prefix_coding.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_qdxt1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_qdxt1.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_qdxt5.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_qdxt5.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_rand.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_rand.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_ray.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_rect.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_resample_filters.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_resample_filters.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_resampler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_resampler.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_rg_etc1.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_rg_etc1.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_ryg_dxt.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_ryg_dxt.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_ryg_types.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_sparse_array.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_sparse_bit_array.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_sparse_bit_array.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_strutils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_strutils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_symbol_codec.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_symbol_codec.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_sysdetection.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_texture_comp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_texture_comp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_texture_conversion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_texture_conversion.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_texture_file_types.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_texture_file_types.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_threaded_clusterizer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_threaded_resampler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_threaded_resampler.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_threading.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_threading_null.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_timer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_timer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_traits.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_tree_clusterizer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_types.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_value.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_value.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_vec.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_vec_interval.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_vector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_vector.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_vector2d.h
    ${CMAKE_CURRENT_SOURCE_DIR}/crn_version.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crnlib.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../inc/crn_decomp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../inc/crnlib.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../inc/dds_defs.h
)

if(WIN32)
    list(APPEND CRNLIB_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/crn_winhdr.h
        ${CMAKE_CURRENT_SOURCE_DIR}/crn_threading_win32.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/crn_threading_win32.h
    )
else()
    list(APPEND CRNLIB_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/crn_threading_pthreads.h
        ${CMAKE_CURRENT_SOURCE_DIR}/crn_threading_pthreads.cpp
    )

    if(APPLE)
        list(APPEND CRNLIB_SRCS
            ${CMAKE_CURRENT_SOURCE_DIR}/crn_darwin_pthreads.h
            ${CMAKE_CURRENT_SOURCE_DIR}/crn_darwin_pthreads.cpp
        )
    endif()
endif()

add_library(crn ${CRNLIB_SRCS})
set_property(TARGET crn PROPERTY CXX_STANDARD 11)
target_include_directories(crn
    PUBLIC
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../inc
    # PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(crn PRIVATE lzma jpge jpgd miniz stb)

include(GenerateExportHeader)
generate_export_header(crn)

if(NOT WIN32)
  find_package(Threads)
  target_link_libraries(crn PUBLIC Threads::Threads)
endif()

install(TARGETS crn)
file(GLOB CRN_INC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../inc/*.h)
list(APPEND CRN_INC_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/crn_export.h)
install(FILES ${CRN_INC_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
