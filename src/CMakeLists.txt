
add_subdirectory(hgp_hrr_code)
add_subdirectory(hgp_vrr_code)


add_library(
    nhfint
    basis.cpp
    eri_class.cpp
    hgp_eri_class.cpp
    eri.cpp
)

target_include_directories(
    nhfint PUBLIC
    .
)

target_link_libraries(
    nhfint PUBLIC
    nhfmath
    hgphrrlib
    hgpvrrlib
)


add_executable(
    main
    main.cpp
)

target_link_libraries(
    main PRIVATE
    nhfint
)