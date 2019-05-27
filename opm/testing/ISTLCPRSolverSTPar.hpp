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
    //typedef Dune::Amg::SequentialInformation POrComm;
    typedef std::size_t GlobalId; // The type for the global index
    typedef Dune::OwnerOverlapCopyCommunication<GlobalId> Communication;
    typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, 1, 1 > > PressureMatrixType;
    typedef Dune::BlockVector< Dune::FieldVector< double, 1 > > PressureVectorType;
    using CoarseOperatorType = Dune::OverlappingSchwarzOperator<PressureMatrixType,
								PressureVectorType,
								PressureVectorType,
								Communication>;
    using FineOperatorType = Dune::OverlappingSchwarzOperator<MatrixType,
							      VectorType,
							      VectorType,
							      Communication>;
    //using FineSmootherType = Dune::SeqILU0<MatrixType, VectorType, VectorType>;
    //using CoarseSmootherType = Dune::SeqILU0<PressureMatrixType, PressureVectorType, PressureVectorType>;
    using FineSmootherType = CprSmootherFine;
    using CoarseSmootherType = CprSmootherCoarse;

    //using Communication =  Dune::Amg::SequentialInformation;
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
      rhs_ = &b;
      //matrix_ = &M;
#if HAVE_MPI
      if( this->isParallel() ) {
	if((parallel_matrix == nullptr) or update_preconditioner){	  
	  parallel_matrix_.reset(new MatrixType(M));
	}else{
	  *parallel_matrix_ = *M; 
	}
	//remove ghost rows in local matrix without doing a copy.
	this->makeOverlapRowsInvalid(*(this->parallel_matrix_));	      
	//Not sure what actual_mat_for_prec is, so put ebosJacIgnoreOverlap as both variables
	//to be certain that correct matrix is used for preconditioning.
	if( ! comm_ or update_preconditioner){	  
	  opAParallel_.reset(new OperatorParallel(*(this->parallel_matrix_), *(this->parallel_matrix_),this->parallelInformation_ ));
	  comm_ = opAParallel_->comm();	  
	  assert(comm_->indexSet().size()==0);
	  //const size_t size = opAParallel_->getmat().N();
	  const size_t size = parallel_matrix_.N();
	  const ParallelISTLInformation& info =
	    boost::any_cast<const ParallelISTLInformation&>( this->parallelInformation_);		    
	  // As we use a dune-istl with block size np the number of components
	  // per parallel is only one.
	  info.copyValuesTo(comm_->indexSet(), comm_->remoteIndices(),
			    size, 1);	  
	}else{
	  opAParallel_.reset(new OperatorParallel(*(this->parallel_matrix_), *(this->parallel_matrix_),comm_));
	}
	using POrComm =  Dune::OwnerOverlapCopyCommunication<int,int>;
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
	constexpr Dune::SolverCategory::Category category=Dune::SolverCategory::overlapping;	      
	auto sp = Dune::createScalarProduct<Vector,POrComm>(*comm_, category);
	sp_ = std::move(sp);
#else
	constexpr int  category = Dune::SolverCategory::overlapping;
	typedef Dune::ScalarProductChooser<Vector, POrComm, category> ScalarProductChooser;
	typedef std::unique_ptr<typename ScalarProductChooser::ScalarProduct> SPPointer;
	SPPointer sp(ScalarProductChooser::construct(info));
	sp_ = std::move(sp);
#endif	
	prepareSolver(*comm_,update_preconditioner);
      }
      else
#endif	 
      {
	if (update_preconditioner) {
	  opAParallel_.reset(new OperatorParallel(*(this->parallel_matrix_), *(this->parallel_matrix_)));
	}
	using POrComm = Dune::Amg::SequentialInformation;
	POrComm parallelInformation_arg;
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
	constexpr Dune::SolverCategory::Category category=Dune::SolverCategory::sequential;
	auto sp = Dune::createScalarProduct<Vector,POrComm>(parallelInformation_arg, category);
	sp_ = std::move(sp);
#else
	constexpr int  category = Dune::SolverCategory::sequential;
	typedef Dune::ScalarProductChooser<Vector, POrComm, category> ScalarProductChooser;
	typedef std::unique_ptr<typename ScalarProductChooser::ScalarProduct> SPPointer;
	SPPointer sp(ScalarProductChooser::construct(parallelInformation_arg));
	sp_ = std::move(sp);
#endif
	prepareSolver(parallelInformation_arg, update_preconditioner);	
      }
    }
    //void prepare(const MatrixType& M, VectorType& b,bool update_preconditioner = true){
    template<typename Comm>  
    void prepareSolver(Comm& comm,bool update_preconditioner = true){
      Vector& istlb = *(this->rhs_);
      comm.copyOwnerToAll(istlb, istlb);      
      if( update_preconditioner or (preconditioner_ == 0) ){
	//.istlMatrix();
	weights_.resize(rhs_->size());
	finesmoother_.reset(new FineSmootherType(*paralel_matrix_,1.0));
	Opm::Amg::getQuasiImpesWeights(*parallel_matrix_, pressureVarIndex, weights_);
	levelTransferPolicy_.reset(new LevelTransferPolicy(comm, weights_));		   
	coarseSolverPolicy_.reset(new CoarseSolverPolicy(smootherArgs_, criterion_, prm_));
	//fineoperator_.reset(new FineOperatorType(*parallel_matrix_));
	preconditioner_.reset(new TwoLevelMethod(opParalel_,//fineoperator_,
						 finesmoother_,
						 *levelTransferPolicy_,
						 *coarseSolverPolicy_,
						 0,
						 1)); 
	linsolve_.reset(new Dune::BiCGSTABSolver<VectorType>(*fineoperator_,
							     *sp_,
							     *preconditioner_,
							     prm_.get<double>("tol"),
							     prm_.get<int>("maxiter"),
							     prm_.get<int>("verbosity")));
	
      }else{
	if(prm_.get<int>("cpr_verbosity")>0){
	  std::cout << " Only update amg solver " << std::endl;
	}
	Opm::Amg::getQuasiImpesWeights(*matrix_, pressureVarIndex, weights_);//set new weights
	finesmoother_.reset(new FineSmootherType(*parallel_matrix_,1.0));//make new finescale smoother
	//levelTransferPolicy_->calculateCoarseEntries(*fineoperator_);// calculate new pressure matrix
	preconditioner_->updatePreconditioner(*fineoperator_,finesmoother_, *coarseSolverPolicy_);//pinfo_);//update coarse solver matrices and smoothers
	//reinterpret_cast<AmgType*>(preconditioner_->updatePreconditioner(opARef, smootherArgs, comm);
      }
    }	  
    
    bool isParallel(){
      return false;
    }
    
    bool solve(VectorType& x) {
      Dune::InverseOperatorResult result;
      VectorType& istlb = *(this->rhs_);
      linsolve_->apply(x, istlb, result);
      checkConvergence(result);
      res_ = result;
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
    std::unique_ptr<MatrixType> parallel_matrix_;// needs a copy to eliminating boundary terms
  //std::unique_ptr< Dune::LinearOperator<Vector, Vector> > opA_;
  
    VectorType weights_;
    pt::ptree prm_;
    typedef typename Dune::Amg::SmootherTraits<CoarseSmootherType>::Arguments  SmootherArgs;
      
    std::shared_ptr<Communication> comm_;
    Criterion criterion_;
    SmootherArgs  smootherArgs_;
    Dune::InverseOperatorResult res_;
    Dune::OverlappingSchwarzScalarProduct<VectorType,Communication> sp_;
    std::shared_ptr<OperatorParallel> opAParallel_;
    //using POrComm =  Dune::OwnerOverlapCopyCommunication<int,int>;
    //std::shared_ptr<POrComm> comm_;
  }; // end ISTLSolver

} // namespace Opm
#endif
