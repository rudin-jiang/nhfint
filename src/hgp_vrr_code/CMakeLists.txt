add_library(
    hgpvrrlib
    boysfun.cpp       hgp_vrr_0_4_.cpp  hgp_vrr_1_4_.cpp  hgp_vrr_2_4_.cpp  hgp_vrr_3_4_.cpp  hgp_vrr_4_4_.cpp
    hgp_vrr_0_0_.cpp  hgp_vrr_1_0_.cpp  hgp_vrr_2_0_.cpp  hgp_vrr_3_0_.cpp  hgp_vrr_4_0_.cpp
    hgp_vrr_0_1_.cpp  hgp_vrr_1_1_.cpp  hgp_vrr_2_1_.cpp  hgp_vrr_3_1_.cpp  hgp_vrr_4_1_.cpp
    hgp_vrr_0_2_.cpp  hgp_vrr_1_2_.cpp  hgp_vrr_2_2_.cpp  hgp_vrr_3_2_.cpp  hgp_vrr_4_2_.cpp
    hgp_vrr_0_3_.cpp  hgp_vrr_1_3_.cpp  hgp_vrr_2_3_.cpp  hgp_vrr_3_3_.cpp  hgp_vrr_4_3_.cpp
)

target_include_directories(
    hgpvrrlib PUBLIC
    .
)