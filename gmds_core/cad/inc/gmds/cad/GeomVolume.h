/*----------------------------------------------------------------------------*/
/** \file    GeomVolume.h
 *  \author  F. LEDOUX
 *  \date    09/21/2010
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_GEOMVOLUME_H_
#define GMDS_GEOM_GEOMVOLUME_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/utils/Exception.h>
#include <gmds/utils/CommonTypes.h>
#include "gmds/cad/GeomEntity.h"

/*----------------------------------------------------------------------------*/
namespace gmds {
    /*----------------------------------------------------------------------------*/
    namespace cad {
        /*----------------------------------------------------------------------------*/
        class GeomPoint;
        class GeomCurve;
        class GeomSurface;
        /*----------------------------------------------------------------------------*/
        /** \class GeomVolume
         *  \brief This class describe the services that are required by the
         *  	   mesh to the geometrical model. As a consequence, this interface only
         *  	   contains query methods.
         */
        /*----------------------------------------------------------------------------*/
        class GeomVolume : public GeomEntity {
        public:
            /*------------------------------------------------------------------------*/
            /** \brief  Constructor
             */
            explicit GeomVolume(const std::string &AName = "Unknown volume")
                : GeomEntity(AName) {
            }

            /*------------------------------------------------------------------------*/
            /** \brief  provides the dimension of the geometrical entity.
             */
            int dim() const override { return 3; }

            /**@brief Check if a point is inside the volume
             */
            virtual bool isIn(const math::Point &AP) const =0;

            /**@brief Accessor to the adjacent points. Warning, there is no
             *  assumption about the ordering
             * @return points that are adjacent to this point
             */
            virtual std::vector<GeomPoint *> &points() =0;

            /**@brief Accessor to the adjacent curves. Warning, there is no
             *  assumption about the ordering
             * @return curves that are adjacent to this point
             */
            virtual std::vector<GeomCurve *> &curves() =0;

            /**@brief Accessor to the adjacent surfaces. Warning, there is no
             *  assumption about the ordering
             * @return surfaces that are adjacent to this point
             */
            virtual std::vector<GeomSurface *> &surfaces() =0;
        };

        /*----------------------------------------------------------------------------*/
    } // namespace cad
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_GEOMVOLUME_H_ */
