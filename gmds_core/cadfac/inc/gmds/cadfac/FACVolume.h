/*----------------------------------------------------------------------------*/
/*
 * FACVolume.h
 *
 *  Created on: 27 juin 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDVOLUME_H_
#define GMDS_GEOM_FACETEDVOLUME_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/cad/GeomVolume.h>
#include <gmds/cadfac/FACSurface.h>
#include <gmds/ig/Mesh.h>

#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> CGAL_Kernel;
typedef CGAL_Kernel::Point_3 CGAL_Point;
typedef CGAL::Surface_mesh<CGAL_Point> CGAL_Surface_Mesh;
/*----------------------------------------------------------------------------*/
#include <tuple>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
/*----------------------------------------------------------------------------*/
/** \class FACVolume
 *  \brief This class implements the volume services that are required by the
 *  	   mesh to the geometrical model.
 *  \param TCoord the data type used to store geometrical data
 */
/*----------------------------------------------------------------------------*/
        class FACVolume : public GeomVolume {

        public:
            /*------------------------------------------------------------------------*/
            /** \brief  Constructor.
            *  @param AMesh mesh support
            *  \param AName	the volume name
            */
            explicit FACVolume(Mesh* AMesh,
                const std::string& AName = "Unknown volume");


            /*------------------------------------------------------------------------*/
            /** \brief  Destructor.
             */
            virtual ~FACVolume();


            /*------------------------------------------------------------------------*/
            /** \brief  computes the area of the entity.
             */
            TCoord computeArea() const override;

            /*------------------------------------------------------------------------*/
            /** \brief  computes the bounding box
             *
             *	\param minXYZ The minimum coordinate of the bounding box.
             *	\param maxXYZ The maximum coordinate of the bounding box.
             */
            void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const override;
				std::tuple<TCoord,TCoord,TCoord,TCoord,TCoord,TCoord>  BBox() const override;

            int id() const override {return m_id;}


            /**@brief Check if a point is inside the volume
             */
            bool isIn(const math::Point& AP) const override;
            /*------------------------------------------------------------------------*/
            /** \brief Project the point AP unto the geometric entity.
            *
            *  \param AP the point to project
             */
            void project(gmds::math::Point& AP) const override;

            /*------------------------------------------------------------------------*/
            /** \brief Get the closest point from AP on the entity
               *  \param AP a 3D point
               *
               *  \return the closest point of APoint on the entity
             */
            math::Point closestPoint(const math::Point& AP) const override;

            /**@brief Accessor to the adjacent points. Warning, there is no
             *  assumption about the ordering
             * @return points that are adjacent to this point
             */
            std::vector<GeomPoint*>& points() override;

            /**@brief Accessor to the adjacent curves. Warning, there is no
             *  assumption about the ordering
             * @return curves that are adjacent to this point
             */
            std::vector<GeomCurve*>& curves() override;

            /**@brief Accessor to the adjacent surfaces. Warning, there is no
             *  assumption about the ordering
             * @return surfaces that are adjacent to this point
             */
            std::vector<GeomSurface*>& surfaces() override;

            /**@brief Reset the global id counter to 1.
             */
             static void resetIdCounter();
        private:

            /** the surface mesh */
            CGAL_Surface_Mesh m_bnd_mesh;
            /** Surface id*/
            int m_id;
            /** global surface index to generate next surface index*/
            static int m_next_id;

            /** adjacent geometric points that bound this volume*/
            std::vector<GeomPoint*> m_adjacent_points;
            /** adjacent geometric curves that bound this volume*/
            std::vector<GeomCurve*> m_adjacent_curves;
            /** adjacent geometric surfaces that bound this volume*/
            std::vector<GeomSurface*> m_adjacent_surfaces;
        };
/*----------------------------------------------------------------------------*/
    } // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_FACETEDVOLUME_H_ */
/*----------------------------------------------------------------------------*/
