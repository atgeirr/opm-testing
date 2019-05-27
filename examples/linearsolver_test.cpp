#include "config.h"

#if !HAVE_ECL_INPUT
#error "The test for the black oil PVT classes requires eclipse input support in opm-common"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
//#include <boost/program_options.hpp>
#include <dune/common/timer.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>        // pt::ptree
#include <boost/property_tree/ini_parser.hpp> 
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
namespace po = boost::program_options;
namespace pt = boost::property_tree;
void printRes(Dune::InverseOperatorResult& res){
  std::cout << "Solver  converged  " << res.converged << std::endl;
  std::cout << "Solver  reduction  " << res.reduction << std::endl;
  std::cout << "Solver  iterations " << res.iterations << std::endl;
  std::cout << "Solver  conv_rate  " << res.conv_rate << std::endl;
}

#include <opm/testing/PressureSolverPolicy.hpp>
#include <opm/testing/PressureTransferPolicy.hpp>
#include <opm/testing/GetQuasiImpesWeights.hpp>
#include <opm/testing/ISTLCPRSolverST.hpp>

void translate_variables_map_to_ptree(po::variables_map &vm, pt::ptree &propTree){

    for(po::variables_map::iterator it=vm.begin(); it!=vm.end(); it++){
        if( it->second.value().type() == typeid(int) ){ propTree.put<int>(it->first,vm[it->first].as<int>()); }
        else if( it->second.value().type() == typeid(float) ){ propTree.put<float>(it->first,vm[it->first].as<float>()); }
        else if( it->second.value().type() == typeid(double) ){ propTree.put<double>(it->first,vm[it->first].as<double>()); }
        else if( it->second.value().type() == typeid(std::string) ){ propTree.put<std::string>(it->first,vm[it->first].as<std::string>()); }
        else if( it->second.value().type() == typeid(size_t) ){ propTree.put<size_t>(it->first,vm[it->first].as<size_t>()); }
        else{ printf("Error: unknown datatype. Abort!\n"); exit(EXIT_FAILURE); }
    }
}



int main(int argc, char **argv)
{
  po::variables_map vm;        
  try {
        po::options_description desc("Allowed options");
        desc.add_options()
	  ("help", "produce help message")
	  ("inputdir", po::value<std::string>(), "input directory")
	  ("outputdir", po::value<std::string>(), "input directory")
	  ("prefix", po::value<std::string>(), "prefix")
	  ("tol", po::value<double>()->default_value(1e-2), "tolerance")
	  ("verbosity", po::value<int>()->default_value(10), "verbosity")
	  ("maxiter", po::value<int>()->default_value(200), "maxiter")
	  ("cpr_verbosity", po::value<int>()->default_value(10), "cpr verbosity")
	  ("cpr_ell_solvetype", po::value<int>()->default_value(2), "cpr solver type")
	  ("cpr_max_ell_iter", po::value<int>()->default_value(1), "cpr maxiter")
        ;

        
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
	  std::cout << desc << "\n";
            return 0;
        }

	if (vm.count("inputdir")) {
	  std::cout << "input directory " 
		    << vm["inputdir"].as<std::string>() << ".\n";
        } else {
	  std::cout << "Input directory not set.\n";
	  return 1;
        }
	if (vm.count("outputdir")) {
	  std::cout << "output directory " 
		    << vm["outputdir"].as<std::string>() << ".\n";
        } else {
	  std::cout << "Output directory not set.\n";
	  return 1;
        }
	if (vm.count("prefix")) {
	  std::cout << "prefix " 
		    << vm["prefix"].as<std::string>() << ".\n";
        } else {
	  std::cout << "Prefix directory not set.\n";
	  return 1;
        }
    }
  catch(std::exception& e) {
      std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
      std::cerr << "Exception of unknown type!\n";
    }
  pt::ptree prm;
  translate_variables_map_to_ptree(vm,prm);
   // if(not(argc == 2)){
  //   std::cout << "Usage xxx dirname" << std::endl;
  //   return 1;
  //}
  std::string dirname(vm["inputdir"].as<std::string>());
  std::string outputdir(vm["outputdir"].as<std::string>());
  std::string inputprefix(dirname);
  inputprefix  += "/";
  inputprefix += vm["prefix"].as<std::string>();
  
  std::string rhsfile(inputprefix + "_rhs_istl.txt");
  std::string matrixfile(inputprefix + "_matrix_istl.txt");
  std::string resultfile(outputdir + "/opm_result");
  Dune::MPIHelper::instance(argc, argv);
  std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
  std::cout << "Rhs file " << rhsfile << std::endl;
  std::cout << "Matrix file " << matrixfile << std::endl;
  typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, 3, 3 > > MatrixType;
  typedef Dune::BlockVector< Dune::FieldVector< double, 3 > > VectorType;
  MatrixType matrix_tmp;
  VectorType rhs_tmp;
  {
    std::ifstream infile(rhsfile);
    if(!infile){
      throw std::runtime_error("Rhs file not read");
    }
    Dune::readMatrixMarket(rhs_tmp,infile);
  }
  {
    std::ifstream infile(matrixfile);
    if(!infile){
      throw std::runtime_error("Matrix file not read");
    }
    Dune::readMatrixMarket(matrix_tmp,infile);
  }
  const MatrixType matrix(matrix_tmp);
  const VectorType rhs(rhs_tmp);
  // std::cout << "Rhs is " << std::endl;
  // Dune::writeMatrixMarket(rhs,std::cout);
  // std::cout << "Matrix is " << std::endl;
  // Dune::writeMatrixMarket(matrix,std::cout);
  std::cout << "Solve system " << std::endl;
  //typedef  SolverType;
  VectorType x(rhs.size());
  {
    Dune::Timer perfTimer;
    perfTimer.start();
    Dune::SeqILU0<MatrixType, VectorType, VectorType> preconditioner(matrix, 1.0);
    Dune::MatrixAdapter<MatrixType, VectorType, VectorType> linearOperator(matrix);
    Dune::BiCGSTABSolver<VectorType> linsolver(linearOperator,
					       preconditioner,
					       vm["tol"].as<double>(), // desired residual reduction factor
					       vm["maxiter"].as<int>(), // maximum number of iterations
					       vm["verbosity"].as<int>()); // verbosity of the solver */
  
  
    //Dune::UMFPack<MatrixType>  linsolver(matrix, 0);
    Dune::InverseOperatorResult res;  
    x = 0.;
    VectorType rhs_loc(rhs);
    linsolver.apply(x, rhs_loc, res);
    double time = perfTimer.stop();
    std::cout << "Solve seqilu0 bicgstab time was " << time << std::endl;
    printRes(res);
    {
      std::ofstream outfile(resultfile + "_bicg_ilu0_" + ".txt");
      if(!outfile){
	throw std::runtime_error("Can not write file");
      }
      Dune::writeMatrixMarket(x,outfile);
    }
  }
  // {
  //   Dune::Timer perfTimer;
  //   perfTimer.start();
  //   Dune::SeqGS<MatrixType, VectorType, VectorType> preconditioner(matrix, 1, 1.0);
  //   Dune::MatrixAdapter<MatrixType, VectorType, VectorType> linearOperator(matrix);
  //   Dune::BiCGSTABSolver<VectorType> linsolver(linearOperator,
  // 					       preconditioner,
  // 					       vm["tol"].as<double>(), // desired residual reduction factor
  // 					       vm["maxiter"].as<int>(), // maximum number of iterations
  // 					       vm["v"].as<int>()); // verbosity of the solver */
  
  
  //   //Dune::UMFPack<MatrixType>  linsolver(matrix, 0);
  //   Dune::InverseOperatorResult res;  
  //   x = 0.;
  //   linsolver.apply(x, rhs, res);
  //   double time = perfTimer.stop();
  //   std::cout << "Solve seqgs bicgstab time was " << time << "  " << std::endl;
  //   printRes(res);
  // }
  {
    //constexpr int pressureEqnIndex = 0;
    constexpr int pressureVarIndex = 1;
    VectorType weights(rhs.size());
    Opm::Amg::getQuasiImpesWeights(matrix, pressureVarIndex, weights);
    Dune::Timer perfTimer;
    perfTimer.start();
    typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, 1, 1 > > PressureMatrixType;
    typedef Dune::BlockVector< Dune::FieldVector< double, 1 > > PressureVectorType;
    using CoarseOperatorType = Dune::MatrixAdapter<PressureMatrixType,PressureVectorType,PressureVectorType>;
    using FineOperatorType = Dune::MatrixAdapter<MatrixType,VectorType,VectorType>;
    using FineSmootherType = Dune::SeqILU0<MatrixType, VectorType, VectorType>;
    using CoarseSmootherType = Dune::SeqILU0<PressureMatrixType, PressureVectorType, PressureVectorType>;
    //using CoarseSmootherType = Dune::SeqGS<PressureMatrixType, PressureVectorType, PressureVectorType>;
    //FineSmootherType finesmoother(matrix, 1.0);
    std::shared_ptr<FineSmootherType> finesmoother = std::make_shared<FineSmootherType>(matrix,1.0);
    
    
    using Criterion  =
      Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<PressureMatrixType,
  								Dune::Amg::FirstDiagonal> >;
    int coarsenTarget=1200;
    int verbose = vm["cpr_verbosity"].as<int>();
    Criterion criterion(15, coarsenTarget);
    criterion.setDebugLevel( verbose ); // no debug information, 1 for printing hierarchy information
    criterion.setDefaultValuesIsotropic(2);
    criterion.setNoPostSmoothSteps( 1 );
    criterion.setNoPreSmoothSteps( 1 );
    //using Smoother = FineSmootherType;
    typedef typename Dune::Amg::SmootherTraits<CoarseSmootherType>::Arguments  SmootherArgs;
    SmootherArgs  smootherArgs;
    smootherArgs.iterations = 1;
    //smootherArgs.relaxationFactor = relax; // seemed to not be set
    // //const Opm::CPRParameter& params(this->parameters_); // strange conversion
    // //ISTLUtility::setILUParameters(smootherArgs, ilu_milu);   
    //smootherArgs.setN(0);
    //smootherArgs.setMilu(0);
    
    using Communication =  Dune::Amg::SequentialInformation;
    Communication comm;
    //std::shared_ptr<Smoother> smoother;
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
    
    LevelTransferPolicy levelTransferPolicy(comm, weights);
    CoarseSolverPolicy coarseSolverPolicy(smootherArgs, criterion, prm);

    FineOperatorType fineoperator(matrix);
    TwoLevelMethod preconditioner(fineoperator,
  				  finesmoother,
  				  levelTransferPolicy,
  				  coarseSolverPolicy,
  				  0,
  				  1);
    


    Dune::MatrixAdapter<MatrixType, VectorType, VectorType> linearOperator(matrix);
    Dune::BiCGSTABSolver<VectorType> linsolver(linearOperator,
  					       preconditioner,
  					       vm["tol"].as<double>(), // desired residual reduction factor
  					       vm["maxiter"].as<int>(), // maximum number of iterations
  					       vm["verbosity"].as<int>()); // verbosity of the solver */
  
  
    //Dune::UMFPack<MatrixType>  linsolver(matrix, 0);
    Dune::InverseOperatorResult res;  
    x = 0.;
    VectorType rhs_loc(rhs);
    linsolver.apply(x, rhs_loc, res);
    double time = perfTimer.stop();
    std::cout << "Solve cpr seqgs bicgstab time was " << time << "  " << std::endl;
    printRes(res);
    {
      std::ofstream outfile(resultfile + "_cpr_ilu0_" + ".txt");
      if(!outfile){
	throw std::runtime_error("Can not write file");
      }
    Dune::writeMatrixMarket(x,outfile);
    }
  }

  {
    Dune::Timer perfTimer;
    perfTimer.start();

    constexpr int pressureVarIndex = 1;
    typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, 1, 1 > > PressureMatrixType;
    typedef Dune::BlockVector< Dune::FieldVector< double, 1 > > PressureVectorType;
    using FineSmootherType = Dune::SeqILU0<MatrixType, VectorType, VectorType>;
    using CoarseSmootherType = Dune::SeqILU0<PressureMatrixType, PressureVectorType, PressureVectorType>;
    Opm::ISTLCprSolverST<MatrixType,
			 VectorType,
			 FineSmootherType,
			 CoarseSmootherType,
			 pressureVarIndex> linsolver(prm);
    VectorType rhs_loc(rhs);
    linsolver.prepare(matrix,rhs_loc);
    x = 0.;
    linsolver.solve(x);
    double time = perfTimer.stop();
    std::cout << "Solve istl cpr seqgs bicgstab time was " << time << "  " << std::endl;
    Dune::InverseOperatorResult res = linsolver.getResult();
    printRes(res);
    {
      std::ofstream outfile(resultfile + "_istl_cpr_ilu0_" + ".txt");
      if(!outfile){
	throw std::runtime_error("Can not write file");
      }
    Dune::writeMatrixMarket(x,outfile);
    }
  }

  //
  std::cout <<"***********************" << std::endl;
  std::cout << "Write result to " << resultfile << std::endl;
  //Dune::writeMatrixMarket(x,std::cout);
  {
    std::ofstream outfile(resultfile + ".txt");
    if(!outfile){
      throw std::runtime_error("Can not write file");
    }
    Dune::writeMatrixMarket(x,outfile);
  }
  std::cout << "Statistics of solve " << std::endl;

  
  // write out option from the program
  // write as jason
  boost::property_tree::json_parser::write_json(std::cout, prm);
}
