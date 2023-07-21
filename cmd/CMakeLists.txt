add_subdirectory(generate_sdf)
add_subdirectory(discrete_field_to_bitmap)
add_subdirectory(generate_density_map)

if (USE_CLIPPER2)
  add_subdirectory(generate_sdf_2d)
  add_subdirectory(discrete_field_to_bitmap_2d)
endif()
