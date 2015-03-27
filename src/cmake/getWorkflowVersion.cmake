
# generate git describe tag and dump to VERSION_FILE
# requires REDIST_DIR and VERSION_FILE

set (GETGIT_CMAKE "${REDIST_DIR}/cmake-modules-c99fd3/GetGitRevisionDescription.cmake")
include ("${GETGIT_CMAKE}")
git_describe(GIT_VERSION --match "v[0-9]*")
if (NOT GIT_VERSION)
    set(GIT_VERSION "UNKNOWN")
endif ()
message(STATUS "Detected workflow version: ${GIT_VERSION}")
file(WRITE ${VERSION_FILE} "${GIT_VERSION}")
