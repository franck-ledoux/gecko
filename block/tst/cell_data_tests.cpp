#include <gecko/block/CellData.h>
#include <gecko/math/BezierCurve.h>
#include <gecko/math/BezierSurface.h>
#include <gecko/math/BezierVolume.h>
#include <gecko/math/Point3d.h>
#include <catch2/catch_test_macros.hpp>

using namespace gecko;

namespace {
    using TestCMap = CMap<BezierCurve<Point3d>, BezierSurface<Point3d>, BezierVolume<Point3d>>;
} // namespace

TEST_CASE("cell_data_cmap_carries_geometry_of_any_degree", "[BlockTestSuite]") {
    // One map type serves every order now: the degree travels in the attributes rather than in the
    // type, which is what lets a structure change order without becoming a different C++ type.
    TestCMap map;
    REQUIRE(map.is_empty());

    const auto edge = map.create_attribute<1>();
    edge->info().curve = BezierCurve<Point3d>(3);
    REQUIRE(edge->info().curve.degree() == 3);
    edge->info().curve = BezierCurve<Point3d>(1);
    REQUIRE(edge->info().curve.degree() == 1);
}
