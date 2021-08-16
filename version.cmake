###############################################################################
# Create configuration header to include version number
###############################################################################

set(GIT_HASH "Unknown")
find_package(Git)
if(GIT_FOUND)
  message(STATUS "Git found")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%h -n 1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )
endif(GIT_FOUND)
configure_file(version_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/version_config.h)
