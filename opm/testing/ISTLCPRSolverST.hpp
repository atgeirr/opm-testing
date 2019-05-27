/*
  Copyright 2016 IRIS AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_ISTLCPRSOLVERST_HEADER_INCLUDED
#define OPM_ISTLCPRSOLVERST_HEADER_INCLUDED

#include <utility>
#include <memory>
#include "PressureTransferPolicy.hpp"
#include "PressureSolverPolicy.hpp"
#include "GetQuasiImpesWeights.hpp"
#include <opm/autodiff/twolevelmethodcpr.hh>
namespace Opm
{
  //=====================================================================
  // Implementation for ISTL-matrix based operator
  //=====================================================================
  
  
  /// This class solves the fully implicit black-oil system by
  /// solving the reduced system (after eliminating well variables)
  /// as a block-structured matrix (one block for all cell variables) for a fixed
  /// number of cell variables np .
  /// \tparam MatrixBlockType The type of the matrix block used.
  /// \tparam VectorBlockType The type of the vector block used.
  /// \tparam pressureIndex The index of the pressure component in the vector
  ///                       vector block. It is used to guide the AMG coarsening.
  ///                       Default is zero.

  template<class MatrixType,
	   class VectorType,
	   class CprSmootherFine,
	   class CprSmootherCoarse,
	   int pressureVarIndex>
  class ISTLCprSolverST  
  {
    typedef typename MatrixType::block_type MatrixBlockType;
    typedef typename VectorType::block_type BlockVector;
    //enum { pressureEqnIndex = BlackOilDefaultIndexTraits::waterCompIdx };
    //enum { pressureVarIndex = Indices::pressureSwitchIdx };
    //typedef WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, false> OperatorSerial;
    typedef Dune::Amg::SequentialInformation POrComm;
    typedef Dune::MatrixAdapter<MatrixType,VectorType, VectorType> MatrixAdapter;
    
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)      
    
#else
    static constexpr int category = Dune::SolverCategory::sequential;
    typedef Dune::ScalarProductChooser<Vector, POrComm, category> ScalarProductChooser;
#endif
    typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, 1, 1 > > PressureMatrixType;
    typedef Dune::BlockVector< Dune::FieldVector< double, 1 > > PressureVectorType;
    using CoarseOperatorType = Dune::MatrixAdapter<PressureMatrixType, PressureVectorType, PressureVectorType>;
    using FineOperatorType = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;
    //using FineSmootherType = Dune::SeqILU0<MatrixType, VectorType, VectorType>;
    //using CoarseSmootherType = Dune::SeqILU0<PressureMatrixType, PressureVectorType, PressureVectorType>;
    using FineSmootherType = CprSmootherFine;
    using CoarseSmootherType = CprSmootherCoarse;
    using Communication =  Dune::Amg::SequentialInformation;
    using Criterion  =
      Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<PressureMatrixType,
								Dune::Amg::FirstDiagonal> >;
    
  public:
    static void registerParameters()
    {
    }

    /// Construct a system solver.
    /// \param[in] parallelInformation In the case of a parallel run
    ///                                with dune-istl the information about the parallelization.
    ISTLCprSolverST(pt::ptree prm)
      : prm_(prm),
	iterations_( 0 ),
	converged_(false),
	criterion_(15,1200)
    {
      criterion_.setDebugLevel(prm_.get<int>("cpr_verbosity")); 
      criterion_.setDefaultValuesIsotropic(2);
      criterion_.setNoPostSmoothSteps( 1 );
      criterion_.setNoPreSmoothSteps( 1 );
      smootherArgs_.iterations = 1; 
      //extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);
      //detail::findOverlapRowsAndColumns(simulator_.vanguard().grid(),overlapRowAndColumns_);	   
    }

    // nothing to clean here
    void eraseMatrix() {
    }

    void prepare(const MatrixType& M, VectorType& b,bool update_preconditioner = true){  
      int newton_iteration = 1;//this->simulator_.model().newtonMethod().numIterations();	  
#if HAVE_MPI			      	  
      if( this->isParallel() )
	{

	}
      else
#endif	    
	{    
	  if( update_preconditioner or (preconditioner_ == 0) ){
	    rhs_ = &b;
	    matrix_ = &M;//.istlMatrix();
	    weights_.resize(rhs_->size());
	    finesmoother_.reset(new FineSmootherType(*matrix_,1.0));
	    Opm::Amg::getQuasiImpesWeights(*matrix_, pressureVarIndex, weights_);
	    levelTransferPolicy_.reset(new LevelTransferPolicy(comm_, weights_));		   
	    coarseSolverPolicy_.reset(new CoarseSolverPolicy(smootherArgs_, criterion_, prm_));
	    fineoperator_.reset(new FineOperatorType(*matrix_));
	    preconditioner_.reset(new TwoLevelMethod(*fineoperator_,
						     finesmoother_,
						     *levelTransferPolicy_,
						     *coarseSolverPolicy_,
						     0,
						     1)); 
	    linsolve_.reset(new Dune::BiCGSTABSolver<VectorType>(*fineoperator_,
								 *preconditioner_,
								 prm_.get<double>("tol"),
								 prm_.get<int>("maxiter"),
								 prm_.get<int>("verbosity")));
		
	  }else{
	    if(prm_.get<int>("cpr_verbosity")>0){
	      std::cout << " Only update amg solver " << std::endl;
	    }
	    Opm::Amg::getQuasiImpesWeights(*matrix_, pressureVarIndex, weights_);//set new weights
	    finesmoother_.reset(new FineSmootherType(*matrix_,1.0));//make new finescale smoother
	    //levelTransferPolicy_->calculateCoarseEntries(*fineoperator_);// calculate new pressure matrix
	    preconditioner_->updatePreconditioner(*fineoperator_,finesmoother_, *coarseSolverPolicy_);//pinfo_);//update coarse solver matrices and smoothers
	  }
	}	  
    }
    bool isParallel(){
      return false;
    }
    
    bool solve(VectorType& x) {
      if( this->isParallel() ){
	// for now only call the superclass
	//bool converged = SuperClass::solve(x);
	//return converged;
      }else{
	// Solve system.
	Dune::InverseOperatorResult result;
	VectorType& istlb = *(this->rhs_);
	linsolve_->apply(x, istlb, result);
	checkConvergence(result);
	res_ = result;
      }
      return this->converged_;
    }
	
    Dune::InverseOperatorResult getResult(){return this->res_;}
    /// Solve the system of linear equations Ax = b, with A being the
    /// combined derivative matrix of the residual and b
    /// being the residual itself.
    /// \param[in] residual   residual object containing A and b.
    /// \return               the solution x

    /// \copydoc NewtonIterationBlackoilInterface::iterations
    int iterations () const { return this->iterations_; }

    /// \copydoc NewtonIterationBlackoilInterface::parallelInformation
    const boost::any& parallelInformation() const { return this->parallelInformation_; }
    void checkConvergence( const Dune::InverseOperatorResult& result ) const
    {
      // store number of iterations
      iterations_ = result.iterations;
      converged_ = result.converged;
	  
      // Check for failure of linear solver.
      if (!result.converged) {
	const std::string msg("Convergence failure for linear solver.");
	//OPM_THROW_NOLOG(NumericalIssue, msg);
      }
    }
    
    void setResidual(VectorType& /* b */) {
      // rhs_ = &b; // Must be handled in prepare() instead.
    }
    
    void getResidual(VectorType& b) const {
      b = *rhs_;
    }
    
    void setMatrix(const MatrixType& /* M */) {
      // matrix_ = &M.istlMatrix(); // Must be handled in prepare() instead.
    }
  protected:      
    // #if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)      
    //       typedef std::shared_ptr< Dune::ScalarProduct<Vector> > SPPointer;
    // #else      
    //       typedef std::unique_ptr<typename ScalarProductChooser::ScalarProduct> SPPointer;
    // #endif		
    std::shared_ptr<FineSmootherType> finesmoother_;
    using LevelTransferPolicy = Opm::PressureTransferPolicy<FineOperatorType,
							    CoarseOperatorType,
							    Communication,
							    pressureVarIndex>;   
    using CoarseSolverPolicy   =
      Opm::Amg::PressureSolverPolicy<CoarseOperatorType,
				     CoarseSmootherType,
				     Criterion,
				     LevelTransferPolicy>;
    using TwoLevelMethod =
      Dune::Amg::TwoLevelMethodCpr<FineOperatorType,
				   CoarseSolverPolicy,
				   FineSmootherType>;

      
    mutable int iterations_;
    mutable bool converged_;
        
      
    std::unique_ptr<LevelTransferPolicy> levelTransferPolicy_;
    std::unique_ptr<CoarseSolverPolicy> coarseSolverPolicy_;
    std::unique_ptr<FineOperatorType> fineoperator_;
    std::unique_ptr<TwoLevelMethod> preconditioner_;
    std::unique_ptr< Dune::BiCGSTABSolver<VectorType> > linsolve_;
    const MatrixType* matrix_;
    VectorType* rhs_;
    VectorType weights_;
    pt::ptree prm_;
    typedef typename Dune::Amg::SmootherTraits<CoarseSmootherType>::Arguments  SmootherArgs;
      
    Communication comm_;
    Criterion criterion_;
    SmootherArgs  smootherArgs_;
    Dune::InverseOperatorResult res_;
  }; // end ISTLSolver

} // namespace Opm
#endif
