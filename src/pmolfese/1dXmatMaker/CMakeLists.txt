# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

cmake_minimum_required(VERSION 3.14)
project(afni_xmat_maker NONE)

# Cargo remains responsible for Rust dependency tracking. An always-invoked
# custom target is cheap when unchanged, and also notices newly added modules.
find_program(CARGO_EXECUTABLE cargo)
if(NOT CARGO_EXECUTABLE)
  message(FATAL_ERROR "1dXmatMaker requires Cargo/Rust 1.99 or newer")
endif()
set(MAKER_TARGET_DIR "${CMAKE_CURRENT_BINARY_DIR}/cargo-target")
set(MAKER_EXECUTABLE "${MAKER_TARGET_DIR}/release/1dXmatMaker")
add_custom_target(1dXmatMaker ALL
  COMMAND "${CARGO_EXECUTABLE}" build
          --manifest-path "${CMAKE_CURRENT_SOURCE_DIR}/Cargo.toml"
          --target-dir "${MAKER_TARGET_DIR}" --release --frozen
  BYPRODUCTS "${MAKER_EXECUTABLE}"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  VERBATIM)

enable_testing()
add_test(NAME maker_library
  COMMAND "${CARGO_EXECUTABLE}" test
          --manifest-path "${CMAKE_CURRENT_SOURCE_DIR}/Cargo.toml"
          --target-dir "${MAKER_TARGET_DIR}" --frozen)
set_tests_properties(maker_library PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
install(PROGRAMS "${MAKER_EXECUTABLE}" DESTINATION bin)

# User documentation and notices accompany the standalone executable.
install(FILES README.md THIRD_PARTY.md DESTINATION share/1dXmatMaker)
install(DIRECTORY licenses/ DESTINATION share/1dXmatMaker/licenses)
install(FILES data/glmsingle/LICENSE DESTINATION share/1dXmatMaker/licenses RENAME GLMsingle.txt)
# The waveform bank is embedded; retain its provenance and redistribution notice.
install(FILES data/glmsingle/LICENSE data/glmsingle/provenance.json
        DESTINATION share/1dXmatMaker/data/glmsingle)
install(DIRECTORY examples DESTINATION share/1dXmatMaker
        FILES_MATCHING PATTERN "*.tsv" PATTERN "*.1D")
