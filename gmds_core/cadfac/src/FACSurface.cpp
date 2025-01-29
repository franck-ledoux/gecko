/*----------------------------------------------------------------------------*/
#include <algorithm>
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
#include "gmds/cadfac/FACSurface.h"
/*----------------------------------------------------------------------------*/
#include "gmds/cadfac/FACCurve.h"
#include "gmds/math/Triangle.h"
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/Mesh.h>

/*----------------------------------------------------------------------------*/
namespace gmds {
	/*----------------------------------------------------------------------------*/
	namespace cad {
		/*----------------------------------------------------------------------------*/
		int FACSurface::m_next_id = 1;
		/*----------------------------------------------------------------------------*/
		void
		FACSurface::resetIdCounter() {
			m_next_id = 1;
		}


		/*----------------------------------------------------------------------------*/
		FACSurface::FACSurface(Mesh *AMesh, std::vector<TCellID> &ADiscret,
		                       const std::string &AName) : GeomSurface(AName), m_support(AMesh), m_id(m_next_id++),
		                                                   m_mesh_faces(ADiscret) {
			buildResearchTree();
		}
		/*----------------------------------------------------------------------------*/
		FACSurface::~FACSurface() {
		}

		/*----------------------------------------------------------------------------*/

		void FACSurface::buildResearchTree() {
			for (auto f_id :m_mesh_faces) {
				auto f = m_support->get<Face>(f_id);
				std::vector<Node> ns = f.get<Node>();
				auto p0 = ns[0].point();
				auto p1 = ns[1].point();
				auto p2 = ns[2].point();
				m_CGAL_triangles.emplace_back(Point(p0.X(), p0.Y(), p0.Z()),
					Point(p1.X(), p1.Y(), p1.Z()),
					Point(p2.X(), p2.Y(), p2.Z()));

			}

			m_AABBTree =  Tree(m_CGAL_triangles.begin(), m_CGAL_triangles.end());
			m_AABBTree.accelerate_distance_queries(); // Optimisation pour les recherches de distances
		}
		/*----------------------------------------------------------------------------*/
		TCoord
		FACSurface::computeArea() const {
			TCoord totalArea = 0.;

			for (unsigned int i = 0; i < m_mesh_faces.size(); i++) {
				totalArea += m_support->get<Face>(m_mesh_faces[i]).area();
			}
			return totalArea;
		}

		/*----------------------------------------------------------------------------*/
		std::tuple<TCoord, TCoord, TCoord, TCoord, TCoord, TCoord>
		FACSurface::BBox() const {
			std::set<TCellID> node_ids;
			for (unsigned int i = 0; i < m_mesh_faces.size(); i++) {
				std::vector<TCellID> f_ids = m_support->get<Face>(m_mesh_faces[i]).getIDs<Node>();
				node_ids.insert(f_ids.begin(), f_ids.end());
			}

			// too many comparisons (3 factor)
			math::Point pi = m_support->get<Node>(*node_ids.begin()).point();
			TCoord minXYZ[3], maxXYZ[3];
			minXYZ[0] = pi.X();
			maxXYZ[0] = pi.X();
			minXYZ[1] = pi.Y();
			maxXYZ[1] = pi.Y();
			minXYZ[2] = pi.Z();
			maxXYZ[2] = pi.Z();
			for (auto i: node_ids) {
				math::Point pi = m_support->get<Node>(i).point();
				minXYZ[0] = std::min(minXYZ[0], pi.X());
				minXYZ[1] = std::min(minXYZ[1], pi.Y());
				minXYZ[2] = std::min(minXYZ[2], pi.Z());
				maxXYZ[0] = std::max(maxXYZ[0], pi.X());
				maxXYZ[1] = std::max(maxXYZ[1], pi.Y());
				maxXYZ[2] = std::max(maxXYZ[2], pi.Z());
			}
			return {minXYZ[0], minXYZ[1], minXYZ[2], maxXYZ[0], maxXYZ[1], maxXYZ[2]};
		}

		/*----------------------------------------------------------------------------*/
		void
		FACSurface::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const {
			auto [min_x, min_y, min_z,max_x,max_y,max_z] = BBox();
			minXYZ[0] = min_x;
			maxXYZ[0] = max_x;
			minXYZ[1] = min_y;
			maxXYZ[1] = max_y;
			minXYZ[2] = min_z;
			maxXYZ[2] = max_z;
		}

		/*----------------------------------------------------------------------------*/

		void
		FACSurface::computeNormal(const math::Point &AP, math::Vector3d &AV) const {
			Face f0 = m_support->get<Face>(m_mesh_faces[0]);
			TCoord min_dist = f0.distance(AP);
			TCellID index = 0;
			for (unsigned int i = 1; i < m_mesh_faces.size(); i++) {
				TCoord dist = m_support->get<Face>(m_mesh_faces[i]).distance(AP);
				if (dist < min_dist) {
					min_dist = dist;
					index = i;
				}
			}
			AV = m_support->get<Face>(m_mesh_faces[index]).normal();
		}

		/*----------------------------------------------------------------------------*/
		math::Point
		FACSurface::closestPoint(const math::Point &AP) const {
			auto closest_point = m_AABBTree.closest_point(Point(AP.X(), AP.Y(), AP.Z()));
			return {closest_point.x(), closest_point.y(), closest_point.z()};
		}

		/*----------------------------------------------------------------------------*/
		math::Point
		FACSurface::center() const {
			throw GMDSException("Not yet implemented");
		}

		/*----------------------------------------------------------------------------*/
		void
		FACSurface::setMeshFaces(const std::vector<Face> &AFaces) {
			m_mesh_faces.clear();
			m_mesh_faces.resize(AFaces.size());
			for (unsigned int i = 0; i < AFaces.size(); i++) {
				m_mesh_faces[i] = AFaces[i].id();
			}
		}

		/*----------------------------------------------------------------------------*/
		void
		FACSurface::getMeshFaces(std::vector<Face> &AFaces) const {
			AFaces.clear();
			AFaces.reserve(m_mesh_faces.size());
			for (auto f_id: m_mesh_faces) {
				AFaces.push_back(m_support->get<Face>(f_id));
			}
		}

		/*----------------------------------------------------------------------------*/
		void
		FACSurface::getTriangulation(std::vector<math::Triangle> &ATri) const {
			ATri.clear();
			ATri.resize(m_mesh_faces.size());
			for (unsigned int i = 0; i < m_mesh_faces.size(); i++) {
				gmds::Face current = m_support->get<gmds::Face>(m_mesh_faces[i]);
				std::vector<Node> nodes = current.get<Node>();
				ATri[i] = math::Triangle(nodes[0].point(), nodes[1].point(), nodes[2].point());
			}
		}

		/*----------------------------------------------------------------------------*/
		void
		FACSurface::propagateOrient(Face AFace, int AMarkTreatedFaces, int AMarkFacesInverted, Mesh *AMesh) {
			std::vector<Face> faces = AFace.get<Face>();

			for (unsigned int iFace = 0; iFace < faces.size(); iFace++) {
				if (!AMesh->isMarked(faces[iFace], AMarkTreatedFaces)) {
					bool faceIsSameOrientation = checkSameOrientFace(AFace, faces[iFace]);

					if (!faceIsSameOrientation) {
						invertFace(faces[iFace]);
						AMesh->mark(faces[iFace], AMarkFacesInverted);
					}

					AMesh->mark(faces[iFace], AMarkTreatedFaces);
					propagateOrient(faces[iFace], AMarkTreatedFaces, AMarkFacesInverted, AMesh);
				}
			}
		}

		/*----------------------------------------------------------------------------*/
		bool
		FACSurface::checkSameOrientFace(Face AFaceRef, Face AFaceCheck) {
			std::vector<Node> nodes = AFaceRef.get<Node>();
			std::vector<Node> nodesBis = AFaceCheck.get<Node>();

			for (unsigned int iNode = 0; iNode < nodes.size(); iNode++) {
				for (unsigned int iNodeBis = 0; iNodeBis < nodesBis.size(); iNodeBis++) {
					if (nodesBis[iNodeBis] == nodes[iNode]) {
						if (nodesBis[(iNodeBis + 1) % nodesBis.size()] == nodes[(iNode + 1) % nodes.size()]) {
							return true;
						} else {
							return false;
						}
					}
				}
			}

			throw GMDSException("FACSurface::checkOrientFace we should not be in this part of the code.");
		}

		/*----------------------------------------------------------------------------*/
		void
		FACSurface::invertAllFaces() {
			std::vector<Face> faces;
			getMeshFaces(faces);
			for (unsigned int iFace = 0; iFace < faces.size(); iFace++) {
				invertFace(faces[iFace]);
			}
		}

		/*----------------------------------------------------------------------------*/
		void
		FACSurface::invertFace(Face AFace) {
			std::vector<Node> nodes = AFace.get<Node>();
			std::vector<Node> nodestmp(nodes.rbegin(), nodes.rend());

			AFace.set<Node>(nodestmp);
		}

		/*----------------------------------------------------------------------------*/
		std::vector<GeomPoint *> &
		FACSurface::points() {
			return m_adjacent_points;
		}

		/*----------------------------------------------------------------------------*/
		std::vector<GeomCurve *> &
		FACSurface::curves() {
			return m_adjacent_curves;
		}

		/*----------------------------------------------------------------------------*/
		std::vector<std::vector<GeomCurve *> >
		FACSurface::loops() {
			std::vector<std::vector<GeomCurve *> > loops;
			// first, we put a boolean flag on each curve to indicate we do not have
			// traversed it
			std::map<int, bool> traversed_curves;
			for (auto c: m_adjacent_curves)
				traversed_curves[c->id()] = false;

			bool remain_curves = true;
			while (remain_curves) {
				remain_curves = false;
				//we pick a non-traversed curve
				int curve_id = -1;
				for (size_t i = 0; i < m_adjacent_curves.size() & remain_curves == false; i++) {
					if (traversed_curves[m_adjacent_curves[i]->id()] == false) {
						remain_curves = true;
						curve_id = i;
					}
				}
				if (remain_curves) {
					//means, we start a new loop
					std::vector<GeomCurve *> cur_loop;
					GeomCurve *seed = m_adjacent_curves[curve_id];
					cur_loop.push_back(seed);
					traversed_curves[seed->id()] = true;
					std::vector<GeomPoint *> front;
					front.insert(front.end(), seed->points().begin(), seed->points().end());
					if (front.size() > 1) {
						// the current loop has two end points and we are going to access to adjacent curves
						// in the boundary surface
						while (!front.empty()) {
							auto pnt = front.back();
							front.pop_back();
							std::vector<GeomCurve *> curves_of_pnt = pnt->curves();
							for (auto c: curves_of_pnt) {
								//we get curves that are adjacent to the current surface only
								std::vector<GeomSurface *> surfaces_of_c = c->surfaces();
								bool find_current_surface = false;
								for (auto s: surfaces_of_c) {
									if (s->id() == this->id())
										find_current_surface = true;
								}
								if (find_current_surface) {
									//this curve bounds the current surface and share a point with the previous curve
									if (!traversed_curves[c->id()]) {
										cur_loop.push_back(c);
										traversed_curves[c->id()] = true;
										std::vector<GeomPoint *> c_pnts = c->points();
										for (auto c_p: c_pnts)
											if (c_p->id() != pnt->id())
												front.push_back(c_p);
									}
								}
							}
						}
					}
					//means the current curve form a loop itself
					loops.push_back(cur_loop);
				}
			}
			return loops;
		}

		/*----------------------------------------------------------------------------*/
		std::vector<GeomVolume *> &
		FACSurface::volumes() {
			return m_adjacent_volumes;
		}

		/*----------------------------------------------------------------------------*/
	} // namespace cad
	/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
