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

    SECTION("Version checks")
    {
        int major = getVersionMajor(env);
        int minor = getVersionMinor(env);
        int patch = getVersionPatch(env);

        CHECK(major == DSP_VERSION_MAJOR);
        CHECK(minor == DSP_VERSION_MINOR);
        CHECK(patch == DSP_VERSION_PATCH);
    }

    SECTION("freeEnv function") {
        freeEnv(env);
        REQUIRE(env == NULL);
    }
}