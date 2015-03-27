
set (THIS_GETGIT_CMAKE "${THIS_REDIST_DIR}/cmake-modules-c99fd3/GetGitRevisionDescription.cmake")
include ("${THIS_GETGIT_CMAKE}")
git_describe(THIS_GIT_VERSION --match "v[0-9]*")
if (NOT THIS_GIT_VERSION)
    set(THIS_GIT_VERSION "UNKNOWN")
endif ()
file(WRITE ${THIS_VERSION_FILE} "${THIS_GIT_VERSION}")
