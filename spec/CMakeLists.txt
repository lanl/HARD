function(spec_headers)
  target_sources(Spec PUBLIC FILE_SET public_headers TYPE HEADERS
    BASE_DIRS ${CMAKE_SOURCE_DIR} FILES ${ARGN})
endfunction()

function(spec_sources)
  target_sources(Spec PRIVATE ${ARGN})
endfunction()

spec_headers(
  control.hh
  io.hh
  limiter.hh
  types.hh
  utils.hh
)
