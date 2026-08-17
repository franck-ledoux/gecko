#include <gecko/block/CellData.h>
#include <gecko/math/BezierCurve.h>
#include <gecko/math/BezierSurface.h>
#include <gecko/math/BezierVolume.h>
#include <gecko/math/Point3d.h>
#include <catch2/catch_test_macros.hpp>

using namespace gecko;

namespace {
    template<std::size_t N>
    using TestCMap = CMap<BezierCurve<N, Point3d>, BezierSurface<N, Point3d>, BezierVolume<N, Point3d>>;
} // namespace

TEST_CASE("cell_data_cmap_instantiates_for_linear_and_curved_order", "[BlockTestSuite]") {
    TestCMap<1> map_linear;
    TestCMap<3> map_cubic;

    REQUIRE(map_linear.is_empty());
    REQUIRE(map_cubic.is_empty());
}
