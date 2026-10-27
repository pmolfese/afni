cmake_minimum_required(VERSION 3.14)
project(1dPetePlot LANGUAGES NONE)
option(PETEPLOT_GUI "Include native windows (macOS/Wayland)" ON)
find_program(CARGO_EXECUTABLE cargo)
if(NOT CARGO_EXECUTABLE)
  message(FATAL_ERROR "1dPetePlot requires Cargo and Rust 1.95 or newer")
endif()
set(PETEPLOT_CARGO_FLAGS --release --locked)
if(NOT PETEPLOT_GUI)
  list(APPEND PETEPLOT_CARGO_FLAGS --no-default-features)
endif()
set(PETEPLOT_CARGO_DIR "${CMAKE_CURRENT_BINARY_DIR}/cargo")
set(PETEPLOT_EXECUTABLE "${PETEPLOT_CARGO_DIR}/release/1dPetePlot${CMAKE_EXECUTABLE_SUFFIX}")
# Run Cargo each time: Cargo tracks source/dependency changes and skips fresh work.
add_custom_target(1dPetePlot ALL
  COMMAND "${CARGO_EXECUTABLE}" build ${PETEPLOT_CARGO_FLAGS}
    --manifest-path "${CMAKE_CURRENT_SOURCE_DIR}/Cargo.toml"
    --target-dir "${PETEPLOT_CARGO_DIR}"
  BYPRODUCTS "${PETEPLOT_EXECUTABLE}"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  VERBATIM USES_TERMINAL)
set(PETEPLOT_INSTALL_BINDIR bin CACHE PATH "Executable installation directory")
install(PROGRAMS "${PETEPLOT_EXECUTABLE}" DESTINATION "${PETEPLOT_INSTALL_BINDIR}")
