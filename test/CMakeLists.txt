add_executable(
    test_basis
    test_basis.cpp
)

target_link_libraries(
    test_basis PRIVATE
    nhfint
    gtest
    gtest_main
)


add_executable(
    test_eri_class
    test_eri_class.cpp
)

target_link_libraries(
    test_eri_class PRIVATE
    nhfint
    gtest
    gtest_main
)