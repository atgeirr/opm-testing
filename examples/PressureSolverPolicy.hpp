#ifndef OPM_PRESSURE_SOLVER_POLICY_HEADER_INCLUDED
#define OPM_PRESSURE_SOLVER_POLICY_HEADER_INCLUDED
#include <dune/istl/paamg/twolevelmethod.hh>
#include <dune/istl/paamg/aggregates.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>
#include <opm/simulators/linalg/amgcpr.hh>
#include <boost/property_tree/ptree.hpp>        // pt::ptree
//#include <boost/property_tree/ini_parser.hpp> 
#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/json_parser.hpp>
//namespace po = boost::program_options;
namespace pt = boost::property_tree;


namespace Opm{
  namespace Amg {
    template<class O, class S, class C, class P>
    class PressureSolverPolicy
    {
    public:
      typedef P LevelTransferPolicy;
      /** @brief The type of the linear operator used. */
      typedef O Operator;
      /** @brief The type of the range and domain of the operator. */
      typedef typename O::range_type X;
      /** @brief The type of the crition used for the aggregation within AMG.*/
      typedef C Criterion;
      /** @brief The type of the communication used for AMG.*/
      typedef typename P::ParallelInformation Communication;
      /** @brief The type of the smoother used in AMG. */
      typedef S Smoother;
      /** @brief The type of the arguments used for constructing the smoother. */
      typedef typename Dune::Amg::SmootherTraits<S>::Arguments SmootherArgs;
      /** @brief The type of the AMG construct on the coarse level.*/
      typedef Dune::Amg::AMGCPR<Operator,X,Smoother,Communication> AMGType;
      /**
       * @brief Constructs the coarse solver policy.
       * @param args The arguments used for constructing the smoother.
       * @param c The crition used for the aggregation within AMG.
       */
      PressureSolverPolicy(const SmootherArgs& args, const Criterion& c, const pt::ptree prm)
	: smootherArgs_(args), criterion_(c), prm_(prm)
      {}
      /** @brief Copy constructor. */
      PressureSolverPolicy(const PressureSolverPolicy& other)
	:  coarseOperator_(other.coarseOperator_), smootherArgs_(other.smootherArgs_),
	   criterion_(other.criterion_), prm_(prm_)
      {}
    private:
      /**
       * @brief A wrapper that makes an inverse operator out of AMG.
       *
       * The operator will use one step of AMG to approximately solve
       * the coarse level system.
       */
      struct AMGInverseOperator : public Dune::InverseOperator<X,X>
      {
        AMGInverseOperator(typename AMGType::Operator& op,
                           const Criterion& crit,
                           const typename AMGType::SmootherArgs& args,
                           const Communication& comm,const pt::ptree& prm)
            :  op_(op), crit_(crit), args_(args), comm_(comm), prm_(prm), amg_()
        {
	  amg_.reset(new AMGType(op, crit, args, comm));
        }

	void updateAmgPreconditioner(typename AMGType::Operator& op){
	  amg_->updateSolver(crit_, op, comm_);
 	}
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
         Dune::SolverCategory::Category category() const override
        {
	  return std::is_same<Communication, Dune::Amg::SequentialInformation>::value ?
	    Dune::SolverCategory::sequential : Dune::SolverCategory::overlapping;
        }
#endif
	

        void apply(X& x, X& b, double reduction, Dune::InverseOperatorResult& res) override
        {
	  DUNE_UNUSED_PARAMETER(reduction);
	  DUNE_UNUSED_PARAMETER(res);
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
	  auto sp = Dune::createScalarProduct<X,Communication>(comm_, op_.category());
#else	  
	  using Chooser = Dune::ScalarProductChooser<X,Communication,AMGType::category>;
	  auto sp = Chooser::construct(comm_);
#endif
	  Dune::Preconditioner<X,X>* prec = amg_.get();
	  // Linear solver parameters
	  const double tolerance = 1e-2;//param_->cpr_solver_tol_;
	  const int maxit        = prm_.get<int>("cpr_max_ell_iter");
	  const int verbosity    = prm_.get<int>("cpr_verbosity");//( param_->cpr_solver_verbose_ &&
				     //comm_.communicator().rank()==0 ) ? 1 : 0;
	  const int solver_type= prm_.get<int>("cpr_ell_solvetype");//
	  if ( solver_type == 0 )
            {
	      // Category of preconditioner will be checked at compile time. Therefore we need
	      // to cast to the derived class
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
	      Dune::BiCGSTABSolver<X> solver(const_cast<typename AMGType::Operator&>(op_), *sp, *prec,
					     tolerance, maxit, verbosity);
#else	      
	      Dune::BiCGSTABSolver<X> solver(const_cast<typename AMGType::Operator&>(op_), *sp,
					     reinterpret_cast<AMGType&>(*prec),
					     tolerance, maxit, verbosity);
#endif
	      solver.apply(x,b,res);
 
            }
	  else if (solver_type == 1)
            {
	      // Category of preconditioner will be checked at compile time. Therefore we need
	      // to cast to the derived class
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
	      Dune::CGSolver<X> solver(const_cast<typename AMGType::Operator&>(op_), *sp, *prec,
				       tolerance, maxit, verbosity);
#else	      
	      Dune::CGSolver<X> solver(const_cast<typename AMGType::Operator&>(op_), *sp,
				       reinterpret_cast<AMGType&>(*prec),
				       tolerance, maxit, verbosity);
#endif
	      solver.apply(x,b,res);
            }
	  else
	    {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
	      Dune::LoopSolver<X> solver(const_cast<typename AMGType::Operator&>(op_), *sp, *prec,
                                         tolerance, maxit, verbosity);
#else	      
	      Dune::LoopSolver<X> solver(const_cast<typename AMGType::Operator&>(op_), *sp,
					 reinterpret_cast<AMGType&>(*prec),
					 tolerance, maxit, verbosity);
#endif
	      solver.apply(x,b,res);	      
	    }
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)

#else
	  delete sp;
#endif
        }

        void apply(X& x, X& b, Dune::InverseOperatorResult& res) override
        {
	  return apply(x,b,1e-8,res);
        }

        ~AMGInverseOperator()
        {}
        AMGInverseOperator(const AMGInverseOperator& other)
	  : x_(other.x_), amg_(other.amg_)
        {
        }
      private:
	//std::unique_ptr<typename AMGType::Operator> op_;
	typename AMGType::Operator& op_;
	Criterion crit_;
	typename AMGType::SmootherArgs args_;
        const Communication& comm_;
	pt::ptree prm_;
        std::unique_ptr<AMGType> amg_;
        X x_;
      };
    public:
      /** @brief The type of solver constructed for the coarse level. */
      typedef AMGInverseOperator CoarseLevelSolver;

      /**
       * @brief Constructs a coarse level solver.
       *
       * @param transferPolicy The policy describing the transfer between levels.
       * @return A pointer to the constructed coarse level solver.
       */
      template<class LTP>
      void setCoarseOperator(LTP& transferPolicy){
	coarseOperator_= transferPolicy.getCoarseLevelOperator();
      }
      template<class LTP>
      CoarseLevelSolver* createCoarseLevelSolver(LTP& transferPolicy)
      {
        coarseOperator_=transferPolicy.getCoarseLevelOperator();
        const LevelTransferPolicy& transfer =
	  reinterpret_cast<const LevelTransferPolicy&>(transferPolicy);
        AMGInverseOperator* inv = new AMGInverseOperator(*coarseOperator_,
                                                         criterion_,
                                                         smootherArgs_,
                                                         transfer.getCoarseLevelCommunication(),
							 prm_);

        return inv; //std::shared_ptr<InverseOperator<X,X> >(inv);

      }
      //void recalculateGalerkin(){
      //	coarseOperator_.recalculateHierarchy();
      //}
    private:
      /** @brief The coarse level operator. */
      std::shared_ptr<Operator> coarseOperator_;
      /** @brief The parameters for the CPR preconditioner. */
      /** @brief The arguments used to construct the smoother. */
      SmootherArgs smootherArgs_;
      /** @brief The coarsening criterion. */
      Criterion criterion_;
      pt::ptree prm_;
      
      //
      
    };
  }
}
#endif
