/*----------------------------------------------------------------------------*/
#include <algorithm>
#include <gmds/cadfac/FACVolume.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace cad {
/*----------------------------------------------------------------------------*/
int FACVolume::m_next_id = 1;
/*----------------------------------------------------------------------------*/
void
FACVolume::resetIdCounter()
{
	m_next_id = 1;
}
/*----------------------------------------------------------------------------*/
FACVolume::FACVolume(Mesh* AMesh,const std::string &AName)
	: GeomVolume(AName),
	  m_id(m_next_id++) {


	std::map<TCellID,CGAL::SM_Vertex_index>  n2n;
	for (auto n_id: AMesh->nodes()) {
		auto n = AMesh->get<Node>(n_id);
		auto p = n.point();
		auto index = m_bnd_mesh.add_vertex({p.X(), p.Y(), p.Z()});
		n2n.insert(std::make_pair(n_id, index));
	}
	//all the faces of bnd_mesh are on the boundary
	for (auto f_id: AMesh->faces()) {
		auto f_node_ids = AMesh->get<Face>(f_id).getIDs<Node>();
		m_bnd_mesh.add_face(n2n[f_node_ids[0]], n2n[f_node_ids[1]],n2n[f_node_ids[2]]);
	}
}

/*----------------------------------------------------------------------------*/
FACVolume::~FACVolume() {}
/*----------------------------------------------------------------------------*/
TCoord
FACVolume::computeArea() const
{
	throw GMDSException("FACVolume::computeArea: not yet implemented");
}
/*----------------------------------------------------------------------------*/
	void FACVolume::project(gmds::math::Point &AP) const {
	throw GMDSException("FACVolume::project: not yet implemented");

	}
	/*----------------------------------------------------------------------------*/
math::Point FACVolume::closestPoint(const math::Point &AP) const {
	throw GMDSException("FACVolume::project: not yet implemented");

}

/*----------------------------------------------------------------------------*/
void
FACVolume::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
{
	auto [min_x, min_y, min_z, max_x, max_y, max_z] = BBox();
	minXYZ[0] = min_x;
	maxXYZ[0] = max_x;
	minXYZ[1] = min_y;
	maxXYZ[1] = max_y;
	minXYZ[2] = min_z;
	maxXYZ[2] = max_z;
}

/*----------------------------------------------------------------------------*/
std::tuple<TCoord, TCoord, TCoord, TCoord, TCoord, TCoord>
FACVolume::BBox() const{
	if (m_adjacent_surfaces.empty()) {
		throw GMDSException("The bounding box of a volume can only be computed if it is enclosed by surfaces!");
	}
	// we have at least one surface
	TCoord minXYZ[3], maxXYZ[3];
	m_adjacent_surfaces[0]->computeBoundingBox(minXYZ, maxXYZ);
	for (size_t s_id = 1; s_id < m_adjacent_surfaces.size(); s_id++) {
		GeomSurface *s = m_adjacent_surfaces[s_id];
		TCoord min_s[3], max_s[3];
		s->computeBoundingBox(min_s, max_s);
		minXYZ[0] = std::min(minXYZ[0], min_s[0]);
		minXYZ[1] = std::min(minXYZ[1], min_s[1]);
		minXYZ[2] = std::min(minXYZ[2], min_s[2]);
		maxXYZ[0] = std::max(maxXYZ[0], max_s[0]);
		maxXYZ[1] = std::max(maxXYZ[1], max_s[1]);
		maxXYZ[2] = std::max(maxXYZ[2], max_s[2]);
	}
	return {minXYZ[0],minXYZ[1],minXYZ[2],maxXYZ[0],maxXYZ[1],maxXYZ[2]};
}
	/*----------------------------------------------------------------------------*/

bool FACVolume::isIn(const math::Point &AP) const {

	CGAL::Side_of_triangle_mesh<CGAL_Surface_Mesh, CGAL_Kernel> inside(m_bnd_mesh);

	return inside({AP.X(),AP.Y(),AP.Z()}) == CGAL::ON_BOUNDED_SIDE;
}
	/*----------------------------------------------------------------------------*/

	std::vector<GeomPoint *> &FACVolume::points()
{
	return m_adjacent_points;
}
/*----------------------------------------------------------------------------*/
std::vector<GeomCurve *> &
FACVolume::curves()
{
	return m_adjacent_curves;
}
/*----------------------------------------------------------------------------*/
std::vector<GeomSurface *> &
FACVolume::surfaces()
{
	return m_adjacent_surfaces;
}
/*----------------------------------------------------------------------------*/
}     // namespace cad
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
