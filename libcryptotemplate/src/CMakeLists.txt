#
# Targets
#

add_library(cryptotemplate 
    PolyCrypto.cpp
    KatePublicParameters.cpp
    Kzg.cpp
    PolyOps.cpp
    Utils.cpp
)

target_include_directories(cryptotemplate PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>"
    "$<INSTALL_INTERFACE:${INSTALL_CMAKE_DIR}/include>"
)

target_link_libraries(cryptotemplate PUBLIC
    ntl ff pthread zm gmp gmpxx xutils xassert
)

include(FindOpenMP)
if(OPENMP_FOUND)
    target_link_libraries(cryptotemplate PUBLIC gomp)
endif()

#target_link_libraries(cryptotemplate PUBLIC Threads::Threads)

#
# Installation
# TODO: Add <Package>Config[Version].cmake files so this package can be easily imported?
# (See https://cmake.org/cmake/help/git-master/manual/cmake-packages.7.html#creating-packages)
#

# This creates the <Package>Config.cmake file and installs it
install(TARGETS cryptotemplate EXPORT cryptotemplateConfig 
   ARCHIVE DESTINATION lib)
install(EXPORT cryptotemplateConfig DESTINATION lib/cmake/cryptotemplate)

# This installs the static or (/and?) dynamic library
install(TARGETS cryptotemplate 
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

# This installs the headers
# WARNING: Don't add / at the end. No slash means cryptotemplate/ directory is created in the destination path
install(DIRECTORY include/cryptotemplate DESTINATION include)
