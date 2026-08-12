#include <gecko/utils/Types.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

using Catch::Approx;
using namespace gecko;
// Définition de types de test
namespace {
    struct TestNodeTag {};
    using TestNodeId = StrongId<TestNodeTag, std::size_t>;

    struct TestFaceTag {};
    using TestFaceId = StrongId<TestFaceTag, std::size_t>;

    using TestSmallId = StrongId<TestNodeTag, std::uint32_t>;
} // namespace

TEST_CASE("StrongId - Initialisation", "[StrongId]") {
    SECTION("Default initialisation ") {
        TestNodeId id;
        REQUIRE_FALSE(id.is_valid());
        REQUIRE(id.value == TestNodeId::invalid_value);
    }

    SECTION("Initialisation with an explicit value") {
        TestNodeId id{42};
        REQUIRE(id.is_valid());
        REQUIRE(id.value == 42);
        REQUIRE(static_cast<std::size_t>(id) == 42);
    }

    SECTION("Reinitialization with reset()") {
        TestNodeId id{100};
        REQUIRE(id.is_valid());
        id.reset();
        REQUIRE_FALSE(id.is_valid());
        REQUIRE(id.value == TestNodeId::invalid_value);
    }

    SECTION("Support of int types (uint32_t)") {
        TestSmallId small_id{123};
        REQUIRE(small_id.is_valid());
        REQUIRE(small_id.value == 123);
        REQUIRE(sizeof(small_id) == sizeof(std::uint32_t));
    }
}

TEST_CASE("StrongId - Comparisons", "[StrongId]") {
    TestNodeId id1{10};
    TestNodeId id2{20};
    TestNodeId id3{10};

    SECTION("Equality and inequality") {
        REQUIRE(id1 == id3);
        REQUIRE(id1 != id2);
    }

    SECTION("Order (<, <=, >, >=)") {
        REQUIRE(id1 < id2);
        REQUIRE(id1 <= id3);
        REQUIRE(id2 > id1);
        REQUIRE(id3 >= id1);
    }
}

TEST_CASE("StrongId - Increment and decrement", "[StrongId]") {
    SECTION("Pre-increment (++id)") {
        TestNodeId id{0};
        TestNodeId res = ++id;

        REQUIRE(id.value == 1);
        REQUIRE(res.value == 1);
    }

    SECTION("Post-increment (id++)") {
        TestNodeId id{0};
        TestNodeId res = id++;

        REQUIRE(id.value == 1);
        REQUIRE(res.value == 0);
    }

    SECTION("Pre-decrement (--id)") {
        TestNodeId id{5};
        TestNodeId res = --id;

        REQUIRE(id.value == 4);
        REQUIRE(res.value == 4);
    }

    SECTION("Post-decrement (id--)") {
        TestNodeId id{5};
        TestNodeId res = id--;

        REQUIRE(id.value == 4);
        REQUIRE(res.value == 5);
    }
}

TEST_CASE("StrongId - Addition and substraction (+, -, +=, -=)", "[StrongId]") {
    TestNodeId id{10};

    SECTION("Operators += et -=") {
        id += 5;
        REQUIRE(id.value == 15);

        id -= 3;
        REQUIRE(id.value == 12);
    }

    SECTION("Operators + et -") {
        TestNodeId res_plus = id + 10;
        REQUIRE(res_plus.value == 20);
        REQUIRE(id.value == 10); // id value must not be updated

        TestNodeId res_minus = id - 4;
        REQUIRE(res_minus.value == 6);
        REQUIRE(id.value == 10);
    }
}
