#include <catch2/catch_test_macros.hpp>
#include <vector>
#include <algorithm>
#include <sstream>
#include <gmds/utils/Exception.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/BitVector.h>

TEST_CASE("UtilsTestSuite - ExceptionWhat") {
    gmds::GMDSException e("test");
    REQUIRE(std::string(e.what()) == "test");
}

TEST_CASE("UtilsTestSuite - common") {
    std::vector<gmds::TCellID> v1 = {1, 2, 3, 4};
    std::vector<gmds::TCellID> v2 = {3, 4, 2, 6};
    std::vector<gmds::TCellID> v3 = gmds::getCommonBut(v1, v2, 4);

    REQUIRE(std::find(v3.begin(), v3.end(), 2) != v3.end());
    REQUIRE(std::find(v3.begin(), v3.end(), 3) != v3.end());
    REQUIRE(std::find(v3.begin(), v3.end(), 4) == v3.end());
}

TEST_CASE("UtilsTestSuite - keepFilter") {
    std::vector<gmds::TCellID> v = {1, 2, 3, 4, 3, 4, 3, 6};
    std::vector<gmds::TCellID> v2 = gmds::keepFilter(v, 2);

    REQUIRE(v2.size() == 2);
    REQUIRE(std::find(v2.begin(), v2.end(), 3) != v2.end());
    REQUIRE(std::find(v2.begin(), v2.end(), 4) != v2.end());

    std::vector<gmds::TCellID> v3 = gmds::keepFilter(v, 3);
    REQUIRE(v3.size() == 1);
    REQUIRE(std::find(v3.begin(), v3.end(), 3) != v3.end());
}

TEST_CASE("UtilsTestSuite - bitVector") {
    gmds::BitVector v;

    REQUIRE(v.size() == 0);
    REQUIRE(v.capacity() == gmds::GChunkSize);
    REQUIRE(v.top() == 0);

    REQUIRE(v.selectNewBit() == 0);
    REQUIRE(v.selectNewBit() == 1);
    REQUIRE(v.selectNewBit() == 2);
    REQUIRE(v.selectNewBit() == 3);
    REQUIRE(v.selectNewBit() == 4);

    int nb_vals = 0;
    for (gmds::TInt val : v) {
        nb_vals++;
    }
    REQUIRE(nb_vals == 5);

    v.unselect(3);
    v.unselect(1);
    REQUIRE(v.size() == 3);

    nb_vals = 0;
    for (auto& val : v) {
        nb_vals++;
    }
    REQUIRE(nb_vals == 3);

    std::stringstream str_in;
    v.serialize(str_in);

    gmds::BitVector v2;
    v2.unserialize(str_in);
    std::stringstream str_out;
    v2.serialize(str_out);

    REQUIRE(str_in.str() == str_out.str());
}

