add_library(
    hgphrrlib
    hgp_hrr_0_0_.cpp
    hgp_hrr_1_0_.cpp
    hgp_hrr_1_1_.cpp
    hgp_hrr_2_0_.cpp
    hgp_hrr_2_1_.cpp
    hgp_hrr_2_2_.cpp
    hgp_hrr_3_0_.cpp
    hgp_hrr_3_1_.cpp
    hgp_hrr_3_2_.cpp
    hgp_hrr_3_3_.cpp
    hgp_hrr_4_0_.cpp
    hgp_hrr_4_1_.cpp
    hgp_hrr_4_2_.cpp
    hgp_hrr_4_3_.cpp
    hgp_hrr_4_4_.cpp
)

target_include_directories(
    hgphrrlib PUBLIC
    .
)