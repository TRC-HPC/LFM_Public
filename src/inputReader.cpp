#include "api/inputReader.h"
#include "api/dictReaderOF.h"
#include "api/fastmesh.h"
#include <dirent.h>
#include <cstring>

bool CInputReader::initialize(int argc, char *argv[])
{
	// Paralell
	m_bParallel = false;

	// Read From Dictionary (system/controlDict)
	CDictReaderOF controlDict("controlDict", "system");
	m_nCommType = controlDict.readDictScalar("commType", "lfm", MPI_EXCHANGE_SPLIT, true);
	m_nHaloCommType = controlDict.readDictScalar("haloCommType", "lfm", MPI_TWOSIDED_NONB, true);
	m_bIsDoublePrecision = controlDict.readDictBool("doublePrecision", "lfm", true, true);
    m_bHaveProbes = controlDict.readDictBool("haveProbes", "lfm/post", false, true);
	m_bHaveSampling = controlDict.readDictBool("haveSampling", "lfm/post", false, true);
	m_bHaveAverage = controlDict.readDictBool("haveAverage", "lfm/post", false, true);
	m_tStartAverage = controlDict.readDictScalar("tStartAverage", "lfm/post", 0, true);
	m_bHaveForces = controlDict.readDictBool("haveForces", "lfm/post", false, true);
	m_bHaveResidual = controlDict.readDictBool("haveResiduals", "lfm/post", false, true);
	m_nSaveForcesSteps = controlDict.readDictScalar("saveForcesStep", "lfm/post", 1, true);
	m_nPrintInfoFrequency = controlDict.readDictScalar("printInfoFreq", "lfm/post", 1, true);
	m_bSaveResiduals = controlDict.readDictBool("saveResiduals", "lfm/post", false, true);
	m_bSaveBlendFactor = controlDict.readDictBool("saveBlendFactor", "lfm/post", false, true);
	m_bSaveRank = controlDict.readDictBool("saveRank", "lfm/post", false, true);
	m_bSaveCellOrder = controlDict.readDictBool("saveCellOrder", "lfm/post", false, true);
	m_dCFLMax = controlDict.readDictScalar("maxCo", "", 1, true);

	// Read From Dictionary (system/fvScheme)
	CDictReaderOF fvScheme("fvSchemes", "system");
	m_nSolver = fvScheme.readDictScalar("solver", "lfm", fastmesh::SOLVER_CAAFOAM, true);
	m_nDimension = fvScheme.readDictScalar("dimension", "lfm", 2, true);
    m_nRkStepCount = fvScheme.readDictScalar("rkOrder", "lfm", 5, true);
	m_bMinmodInterpolation = fvScheme.readDictBool("minmodExists", "lfm", false, true);
	m_bConstantTimeStep = fvScheme.readDictBool("constantTimeStep", "lfm", true, true);

	// Read From Dictionary (constant/spongeDict)
    CDictReaderOF spongeDict("spongeDict", "constant");			
	m_dPInf = spongeDict.readDictScalar("pinf");
	m_dTInf = spongeDict.readDictScalar("Tinf");
	m_dUInf[0] = spongeDict.readDictScalar("Uinf_x", "", 0, true);
	m_dUInf[1] = spongeDict.readDictScalar("Uinf_y", "", 0, true);
	m_dUInf[2] = spongeDict.readDictScalar("Uinf_z", "", 0, true);
	m_dLs = spongeDict.readDictScalar("Ls");
	m_dMach = spongeDict.readDictScalar("M");
	m_dK = spongeDict.readDictScalar("k", "", 4.34294481903252e-01, true);

	// Read From Dictionary (constant/thermophysicalProperties)
	CDictReaderOF thermophysicalProperties("thermophysicalProperties", "constant");
	m_dCp = thermophysicalProperties.readDictScalar("Cp", "mixture/thermodynamics");
	m_dMolWeight = thermophysicalProperties.readDictScalar("molWeight", "mixture/specie");
	m_dMu0 = thermophysicalProperties.readDictScalar("mu", "mixture/transport");
	m_dPrandtl = thermophysicalProperties.readDictScalar("Pr", "mixture/transport");

	// Read From Dictionary (constant/turbulenceProperties)
	CDictReaderOF turbulenceProperties("turbulenceProperties","constant");
	std::string m_sSim_type = turbulenceProperties.readDictString("simulationType");
	if(m_sSim_type.compare("laminar")==0){
		m_bisLaminar = true;
	} else{
		m_bisLaminar = false;
	}

	// Read RunTime Arguments
	
	return (readRuntimeArgs(argc, argv) == 1);
}

// -------------------------------------------------------------------------- //
void CInputReader::showUsage( char *argv[] )
{
    std::cerr << "Usage: " << argv[0] << " <option(s)> ARGUMENT\n"
              << "Options:\n"
              << "\t-h,--help          \tShow this help message\n"
			  << "\t-parallel		    \tRun in parallel\n"
			  << "\t-c,--comm_type     \tcommunicator type              \tSpecify MPI communicator type\n"
			  << "\t\t available communicators: onesided_nonb, onesided_blck, twosided_nonb, twosided_blck, neighcoll_nonb, neighcoll_blck, neighcoll_pers"
			  << "\tLFM Expected Data:\n"
			  << "\t\t-polyMesh grids: located under processor#rank/constant/polyMesh\n"			  
			  << "\t\t-system/controlDict\n"
			  << "\t\t-system/fvScheme\n"
			  << "\t\t-constant/spongeDict\n"
			  << "\t\t-constant/thermophysicalProperties\n"			  
              << std::endl;
}

// -------------------------------------------------------------------------- //
bool CInputReader::parseHaloCommType(const char* sCommType)
{
	if(strcmp(sCommType, "onesided_nonb") == 0)
		m_nHaloCommType = MPI_ONESIDED_NONB;
	else if(strcmp(sCommType, "onesided_blck") == 0)
		m_nHaloCommType = MPI_ONESIDED_BLCK;
	else if(strcmp(sCommType, "twosided_nonb") == 0)
		m_nHaloCommType = MPI_TWOSIDED_NONB;
	else if(strcmp(sCommType, "twosided_blck") == 0)
		m_nHaloCommType = MPI_TWOSIDED_BLCK;
	else if (strcmp(sCommType, "twosided_pers") == 0)
		m_nHaloCommType = MPI_TWOSIDED_PERS;
	else if(strcmp(sCommType, "neighcoll_nonb") == 0)
		m_nHaloCommType = MPI_NEIGHCOLL_NONB;
	else if(strcmp(sCommType, "neighcoll_blck") == 0)
		m_nHaloCommType = MPI_NEIGHCOLL_BLCK;
	else if(strcmp(sCommType, "neighcoll_pers") == 0)
		m_nHaloCommType = MPI_NEIGHCOLL_PERS;
	else
		return false;
	return true;
}

// -------------------------------------------------------------------------- //
bool CInputReader::readRuntimeArgs(int argc, char *argv[])
{
	// No Arguments
	if (argc < 2)
		return true;

	int nRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
	const bool bIsMaster = (nRank == 0);
		
	// Go over options
	while (true) 
	{		
		// Define required arguments
		int option_index = 0;
		static struct option long_options[] = {
			{"help", no_argument, 0, 0 },
			{"parallel", no_argument, 0, 0},
			{"comm_type", required_argument, 0, 0 },
		};

		// Read argument
		int c = getopt_long(argc, argv, "hpc:", long_options, &option_index);

		if (c == -1) break;

		switch (c) {
			// Long option
			case 0:
            	if(option_index == 0 ){
					showUsage(argv);						
					return false;
				}
				else if(option_index == 1)
					m_bParallel = true;					
				else if( option_index == 2 ){
					if(!parseHaloCommType(optarg))
					{
						if(bIsMaster)
						{
							std::cerr << "unknown communicator type: " << optarg << std::endl;
							showUsage(argv);						
						}
						return false;
					}
				}
            	break;

			// Short options (if needed)
			case 'c':
				if(!parseHaloCommType(optarg))
				{
					if(bIsMaster)
					{
						std::cerr << "unknown communicator type: " << optarg << std::endl;
						showUsage(argv);						
					}
					return false;
				}
				break;
			case 'p':
				m_bParallel = true;
				break;
			default:
				if(bIsMaster)
					showUsage(argv);						
				return false;
        }
	}

	// If received invalid runTime argument
	if (optind < argc) {
		if(bIsMaster)
		{
			printf("non-option ARGV-elements: ");
			while (optind < argc)
				printf("%s ", argv[optind++]);
			printf("\n");
		}
	}

	return true;
}
