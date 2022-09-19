#include <gtest/gtest.h>
#include <iostream>
#include <random>
#include "eri_class.hpp"

static const double absErr = 10e-12;

TEST(TestEriClass, TestConstructorWithoutEri) {
    for (std::size_t a = 0; a <= 6; ++a) {
    for (std::size_t b = 0; b <= 6; ++b) {
    for (std::size_t c = 0; c <= 6; ++c) {
    for (std::size_t d = 0; d <= 6; ++d) {
        nhfInt::EriClass ec(a, b, c, d);
        
        EXPECT_EQ(ec.angA, a);
        EXPECT_EQ(ec.angB, b);
        EXPECT_EQ(ec.angC, c);
        EXPECT_EQ(ec.angD, d);

        EXPECT_EQ(ec.nbsA, (a+1)*(a+2)/2);
        EXPECT_EQ(ec.nbsB, (b+1)*(b+2)/2);
        EXPECT_EQ(ec.nbsC, (c+1)*(c+2)/2);
        EXPECT_EQ(ec.nbsD, (d+1)*(d+2)/2);

        std::size_t size = ec.nbsA * ec.nbsB 
                         * ec.nbsC * ec.nbsD;
        EXPECT_EQ(ec.eri.size(), size);

        for (std::size_t i = 0; i < ec.eri.size(); ++i){
            EXPECT_EQ(ec.eri[i], 0);
        }

        std::size_t num2 = ec.nbsB * ec.nbsC * ec.nbsD;
        std::size_t num3 = ec.nbsC * ec.nbsD;
        std::size_t num4 = ec.nbsD;

        EXPECT_EQ(ec.num2, num2);
        EXPECT_EQ(ec.num3, num3);
        EXPECT_EQ(ec.num4, num4);
    }}}}
}

TEST(TestEriClass, TestFourIndxAccessOperator) {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    for (std::size_t a = 0; a <= 6; ++a) {
    for (std::size_t b = 0; b <= 6; ++b) {
    for (std::size_t c = 0; c <= 6; ++c) {
    for (std::size_t d = 0; d <= 6; ++d) {
        
        std::size_t nbsA = (a+1) * (a+2) / 2;
        std::size_t nbsB = (b+1) * (b+2) / 2;
        std::size_t nbsC = (c+1) * (c+2) / 2;
        std::size_t nbsD = (d+1) * (d+2) / 2;

        std::size_t size = nbsA * nbsB * nbsC * nbsD;
        std::vector<double> eri(size);

        for (std::size_t i = 0; i < eri.size(); ++i) {
            eri[i] = dis(gen);
        }

        nhfInt::EriClass ec(a, b, c, d, eri);

        std::size_t pos = 0;
        for (std::size_t i = 0; i < ec.nbsA; ++i) {
        for (std::size_t j = 0; j < ec.nbsB; ++j) {
        for (std::size_t k = 0; k < ec.nbsC; ++k) {
        for (std::size_t l = 0; l < ec.nbsD; ++l) {
            EXPECT_NEAR(ec(i, j, k, l), eri[pos++], absErr);
        }}}}
    }}}}
}

TEST(TestEriClass, TestTwoIndxAccessOperator) {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    for (std::size_t a = 0; a <= 6; ++a) {
    for (std::size_t b = 0; b <= 6; ++b) {
    for (std::size_t c = 0; c <= 6; ++c) {
    for (std::size_t d = 0; d <= 6; ++d) {
        
        std::size_t nbsA = (a+1) * (a+2) / 2;
        std::size_t nbsB = (b+1) * (b+2) / 2;
        std::size_t nbsC = (c+1) * (c+2) / 2;
        std::size_t nbsD = (d+1) * (d+2) / 2;

        std::size_t size = nbsA * nbsB * nbsC * nbsD;
        std::vector<double> eri(size);

        for (std::size_t i = 0; i < eri.size(); ++i) {
            eri[i] = dis(gen);
        }

        nhfInt::EriClass ec(a, b, c, d, eri);

        std::size_t pos = 0;
        for (std::size_t i = 0; i < ec.nbsA * ec.nbsB; ++i) {
        for (std::size_t j = 0; j < ec.nbsC * ec.nbsD; ++j) {
            EXPECT_NEAR(ec(i, j), eri[pos++], absErr);
        }}
    }}}}
}
