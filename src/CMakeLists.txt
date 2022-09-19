add_library(
    nhfint
    basis.cpp
    eri_class.cpp
)

target_include_directories(
    nhfint PUBLIC
    .
)

target_link_libraries(
    nhfint PUBLIC
    nhfmath
)


add_executable(
    main
    main.cpp
)

target_link_libraries(
    main PRIVATE
    nhfint
)