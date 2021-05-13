set(cryptotemplate_bench_sources
    BenchConvertAndMultiexp.cpp
    BenchExp.cpp
    BenchFk.cpp
    BenchFFT.cpp
    BenchMultiexp.cpp
    BenchNtlConversion.cpp
    BenchPairing.cpp
    BenchPointproofs.cpp
    BenchPolynomialOps.cpp
)

foreach(appSrc ${cryptotemplate_bench_sources})
    get_filename_component(appName ${appSrc} NAME_WE)
    set(appDir ../bin/bench)

    add_executable(${appName} ${appSrc})
    target_link_libraries(${appName} PRIVATE cryptotemplate)

    set_target_properties(${appName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${appDir})
endforeach()
