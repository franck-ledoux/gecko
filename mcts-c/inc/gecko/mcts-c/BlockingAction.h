/*----------------------------------------------------------------------------*/
#ifndef GMDS_MCTS_BLOCKING_ACTION_H
#define GMDS_MCTS_BLOCKING_ACTION_H
/*----------------------------------------------------------------------------*/
#include "mcts/IAction.h"
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gecko/cblock/Blocking.h>
/*----------------------------------------------------------------------------*/
namespace gecko {
/*----------------------------------------------------------------------------*/
namespace mctsc {
/*----------------------------------------------------------------------------*/
/**@brief this class encapsulate the action of cutting a block edge
 */
class EdgeCutAction : public IAction {
 public:
	/**@brief Computes the state reach from @p AState by applying the current
     * action
     * @param[in] AState the state we start from
     * @return the state that is built from this @p AState when applying the
     * action
	 */
	std::shared_ptr<IState> apply_on(std::shared_ptr<IState> AState) const override;

	EdgeCutAction(const gmds::TCellID AEdgeId, const double AParam,const gmds::math::Point APoint);

	bool operator==(const IAction& other) const override;
	std::string get_description() const override;

 private:
	/** the id of the edge to cut */
	gmds::TCellID m_edge_id;
	/** the parameter where the edge is cutted */
	double m_cut_param;
	/** the point that we try to capture */
	gmds::math::Point m_capt_point;
};

/*----------------------------------------------------------------------------*/
/**@brief this class encapsulate the action of cutting a block edge
 */
class BlockRemovalAction : public IAction {
 public:
	/**@brief Computes the state reach from @p AState by applying the current
     * action
     * @param[in] AState the state we start from
     * @return the state that is built from this @p AState when applying the
     * action
	 */
	std::shared_ptr<IState> apply_on(std::shared_ptr<IState> AState) const override;

	BlockRemovalAction(const gmds::TCellID ABlockID);

	bool operator==(const IAction& other) const override;
	std::string get_description() const override;

 private:
	/** the id of the block to remove*/
	gmds::TCellID m_block_id;
};

/*----------------------------------------------------------------------------*/
/**@brief this class encapsulate the action of trying to classify boundary edges
 * and faces
 */
class CaptureAction : public IAction {
 public:
	/**@brief Computes the state reach from @p AState by applying the current
     * action
     * @param[in] AState the state we start from
     * @return the state that is built from this @p AState when applying the
     * action
	 */
	std::shared_ptr<IState> apply_on(std::shared_ptr<IState> AState) const override;

	CaptureAction();

	bool operator==(const IAction& other) const override;
	std::string get_description() const override;

};
/*----------------------------------------------------------------------------*/
}     // namespace blocking
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MCTS_BLOCKING_ACTION_H
