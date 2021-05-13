set(cryptotemplate_test_sources
    TestGroupFft.cpp
    TestKatePublicParams.cpp
    TestPolyOps.cpp
    TestLibff.cpp
    TestParallelPairing.cpp
    TestRootsOfUnity.cpp
    TestXGCD.cpp
)

foreach(appSrc ${cryptotemplate_test_sources})
    get_filename_component(appName ${appSrc} NAME_WE)
    set(appDir ../bin/test)

    add_executable(${appName} ${appSrc})
    target_link_libraries(${appName} PRIVATE cryptotemplate)

    add_test(NAME ${appName} COMMAND ${appName})
    set_target_properties(${appName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${appDir})
endforeach()
