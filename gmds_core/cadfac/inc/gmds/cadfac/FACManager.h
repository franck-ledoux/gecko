/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDGEOMMANAGER_H_
#define GMDS_GEOM_FACETEDGEOMMANAGER_H_
/*----------------------------------------------------------------------------*/
#include <stack>
/*----------------------------------------------------------------------------*/
// gmds File headers
/*----------------------------------------------------------------------------*/
#include "gmds/utils/CommonTypes.h"
#include "gmds/utils/Exception.h"

#include "gmds/math/Vector.h"

#include "gmds/cad/GeomCurve.h"
#include "gmds/cad/GeomManager.h"
#include "gmds/cadfac/FACCurve.h"
#include "gmds/cadfac/FACPoint.h"
#include "gmds/cadfac/FACSurface.h"
#include "gmds/cadfac/FACVolume.h"

#include "gmds/cad/GeomMeshLinker.h"
#include "gmds/ig/Mesh.h"
#include "gmds/ig/MeshDoctor.h"

/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace cad {
/*----------------------------------------------------------------------------*/
/** \class  FACManager
 *  \brief  This class creates all the model entities and services relative to
 *  		the faceted representation of the geometry.
 *
 *
 */
/*----------------------------------------------------------------------------*/
class FACManager : public GeomManager
{

 public:
	/*------------------------------------------------------------------------*/
	/** \brief  Default constructor
	 */
	FACManager();

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor
	 */
	virtual ~FACManager();

	/*------------------------------------------------------------------------*/
	/** @brief create a faceted model from a 3D tet mesh. There is no assumption
	 *         about AFromMesh, we get its boundary to build the model.
	 *
	 *         Warning, it means that we totally rebuild the model and we
	 *         consider that it is made of only one single volume!
	 *
	 *  @param AFromMesh 3D mesh used to build the model
	 */
	void initFrom3DMesh(Mesh *AFromMesh);

	/*------------------------------------------------------------------------*/
	/** @brief create a faceted model from a 3D tet mesh. There is no assumption
	 *         about AFromMesh, we get its boundary to build the model. A linker
	 *         given in arguments is filled up to keep the classification.
	 *
	 *         Warning, it means that we totally rebuild the model and we
	 *         consider that it is made of only one single volume!
	 *
	 *  @param AFromMesh 3D mesh used to build the model
	 */
	void initAndLinkFrom3DMesh(Mesh *AFromMesh, GeomMeshLinker *ALinker);

	/*------------------------------------------------------------------------*/
	/** \brief  reinitializes the complete model from data stored in mesh_.
	 */
	void updateFromMesh();

	/*------------------------------------------------------------------------*/
	/** @brief reinitialize the geom model
	 */
	void clear();

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of points of the model.
	 *
	 *	\return the number of points.
	 */
	 TInt getNbPoints() const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of curves of the model.
	 *
	 *	\return the number of curves.
	 */
	TInt getNbCurves() const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of surfaces of the model.
	 *
	 *	\return the number of surfaces.
	 */
	TInt getNbSurfaces() const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of volumes of the model.
	 *
	 *	\return the number of volumes.
	 */
	TInt getNbVolumes() const override;

	GeomEntity *getEntity(TInt AID, TInt ADim) override;
	/*------------------------------------------------------------------------*/
	/** \brief  Access to the points of the model.
	 *
	 *  \param points the points of the model.
	 */
	void getPoints(std::vector<GeomPoint *> &points) const override;
	std::vector<GeomPoint*> getPoints() const override;
	/*------------------------------------------------------------------------*/
	/** \brief  Access to the curves of the model.
	 *
	 *  \param curves the curves of the model.
	 */
	void getCurves(std::vector<GeomCurve *> &curves) const override;
	std::vector<GeomCurve *> getCurves() const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the surface of the model.
	 *
	 *  \param surfaces the surfaces of the model.
	 */
	void getSurfaces(std::vector<GeomSurface *> &surfaces) const override;
	std::vector<GeomSurface *> getSurfaces() const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the volumes of the model.
	 *
	 *  \param volumes the volumes of the model.
	 */
	void getVolumes(std::vector<GeomVolume *> &volumes) const override;
	std::vector<GeomVolume *> getVolumes() const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Gives access to the point of id @AID, Return NullPtr if it does
	 *          not exist.
	 *  \return A point
	 */
	GeomPoint *getPoint(TInt AID) override;
	/*------------------------------------------------------------------------*/
	/** \brief  Gives access to the curve of id @AID, Return NullPtr if it does
	 *          not exist.
	 *  \return A curve
	 */
	 GeomCurve *getCurve(TInt AID) override;
	/*------------------------------------------------------------------------*/
	/** \brief  Gives access to the surface of id @AID, Return NullPtr if it
	 *          does not exist.
	 *  \return A surface
	 */
	 GeomSurface *getSurface(TInt AID) override;
	/*------------------------------------------------------------------------*/
	/** \brief  Gives access to the volume of id @AID, Return NullPtr if it
	 *          does not exist.
	 *  \return A volume
	 */
	GeomVolume *getVolume(TInt AID) override;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the curve common to 2 points
	 *
	 *  \param return the id of the common curve, and -1 if it doesn't exist
	 */
	int getCommonCurve(GeomPoint *AP1, GeomPoint *AP2) const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the surfaces common to 2 points
	 *  \param AP1 first point
	 * \param AP2 second point
	 *
	 *  \return the ids of the common surfaces (potentially empty)
	 */
	std::vector<int> getCommonSurfaces(GeomPoint *AP1, GeomPoint *AP2) const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the surface common to 2 curves
	 *  \param AC1 first curve
	 * \param AC2 second curve
	 *  \return the id of the common surface, and -1 if it doesn't exist
	 */
	int getCommonSurface(GeomCurve *AC1, GeomCurve *AC2) const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Write the mesh representation
	 */
	void write_surfaces(std::string AFilename) const;


 private:
	/*------------------------------------------------------------------------*/
	/** \brief  Build topological relations between geom entities. Must be done
	 *          in the init functions only!
	 */
	void buildTopologicalConnectivity(const int dim = 3);

	void OutwardNormal(Face *f);
	bool IsInsideFace(const int &, const math::Vector3d &, const math::Vector3d &, Face fi);

	math::Vector3d GetBarycentricCoefs(math::Vector3d &V1, math::Vector3d &V2, math::Vector3d &V3, math::Vector3d &I, math::Vector3d &normal);

	Mesh m_mesh;

	std::vector<FACVolume *> m_volumes;
	std::vector<FACSurface *> m_surfaces;
	std::vector<FACCurve *> m_curves;
	std::vector<FACPoint *> m_points;

	std::map<TInt, TInt> m_map_node_var_to_point;
	std::map<TInt, TInt> m_map_edge_var_to_curve;
	std::map<TInt, TInt> m_map_face_var_to_surf;
};
/*----------------------------------------------------------------------------*/
}     // namespace cad
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* gmds_GEOM_FACETEDGEOMMANAGER_H_ */
/*----------------------------------------------------------------------------*/
