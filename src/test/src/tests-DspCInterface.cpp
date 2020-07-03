// tests-BlkModel.cpp
#include "catch.hpp"

#include "DspCInterface.h"

TEST_CASE("C interface functions") {
    DspApiEnv* env = NULL;

    SECTION("createEnv function") {
        env = createEnv();
        REQUIRE(env != NULL);
        REQUIRE(env->model_ == NULL);
        REQUIRE(env->solver_ == NULL);
    }

    SECTION("freeEnv function") {
        freeEnv(env);
        REQUIRE(env == NULL);
    }
}