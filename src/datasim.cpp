/*****************************************************************************
*    <DataSim: VLBI data simulator>                                          *
*    Copyright (C) <2015> <Zheng Meyer-Zhao>                                 *
*                  <2026> <Aard Keimpema>
*                                                                            *
*    This file is part of DataSim.                                           *
*                                                                            *
*    DataSim is free software: you can redistribute it and/or modify         *
*    it under the terms of the GNU General Public License as published by    *
*    the Free Software Foundation, either version 3 of the License, or       *
*    (at your option) any later version.                                     *
*                                                                            *
*    DataSim is distributed in the hope that it will be useful,              *
*    but WITHOUT ANY WARRANTY; without even the implied warranty of          *
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
*    GNU General Public License for more details.                            *
*                                                                            *
*    You should have received a copy of the GNU General Public License       *
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
*****************************************************************************/

#include <iostream>
#include <cstddef>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include <getopt.h>
#include <mpi.h>
#include "datasim.h"
#include "jsonconfig.h"
#include "delaytable.h"
#include "util.h"
#include "subband.h"
#include "vdifio.h"
#include "vdifzipper.h"
#include "catvdif.h"

using namespace std;

static void usage(int argc, char **argv)
{
  cout << endl;
  cout << "Usage:  " << argv[0] << " [<options>] <config.json>" << endl;
  cout << endl;
  cout << "  All simulation settings (start time, duration, stations, SEFDs,\n"
       << "  SFXC delay tables, bands, ...) are read from the JSON configuration\n"
       << "  file; see the examples directory for the schema." << endl;
  cout << endl;
  cout << "  options can include:" << endl;
  cout << "     -h" << endl;
  cout << "     --help        display this information and quit." << endl;
  cout << endl;
  cout << "     -v" << endl;
  cout << "     --verbose     increase the verbosity of the output; -v -v for more." << endl;
  cout << endl;
  cout << "     -t" << endl;
  cout << "     --test        run in test mode, generate 1 second data for each station,\n"
       << "                   no matter what's given in the configuration file." << endl;
  cout << endl;
}

static void cmdparser(int argc, char* argv[], Setup &setupinfo)
{
  char tmp;
  static struct option long_options[] = {
    {"help",      no_argument,        0,  'h'},
    {"verbose",   no_argument,        0,  'v'},
    {"test",      no_argument,        0,  't'},
    {0,           0,                  0,   0 }
  };
  int long_index = 0;
  while((tmp=getopt_long(argc,argv,"hvt",
              long_options, &long_index)) != -1)
  {
    switch(tmp)
    {
      case 'h':
        usage(argc, argv);
        exit (EXIT_SUCCESS);
        break;
      case 'v':
        setupinfo.verbose++;
        break;
      case 't':
        setupinfo.test = 1;
        break;
      default:
        usage(argc, argv);
        exit (EXIT_FAILURE);
        break;
    }
  }
  // check if a JSON configuration file name is given
  if(optind > argc - 1)
  {
    cerr << "No JSON configuration file provided, nothing to do ..." << endl;
    usage(argc, argv);
    exit (EXIT_FAILURE);
  }
  else if(optind < argc - 1)
  {
    cerr << "Multiple input files provided, only one is expected." << endl;
    usage(argc, argv);
    exit (EXIT_FAILURE);
  }
  else
  {
    strcpy(setupinfo.inputfilename, argv[optind]);
  }
}

/*
 * Fill the simulation settings of the setup struct from the JSON
 * configuration (every rank parses the JSON file itself).
 */
static void applyconfig(JsonConfig* config, Setup &setupinfo)
{
  int numdatastreams = config->getNumDataStreams();
  if(numdatastreams > MAXANT)
    throw runtime_error("more stations in JSON configuration than supported (MAXANT)");

  setupinfo.seed = config->getSeed();
  setupinfo.sfluxdensity = config->getFluxDensity();
  setupinfo.numdivs = config->getNumDivs();
  setupinfo.pcal = config->getPcal();
  setupinfo.specres = config->getSpecResScale();
  for(int i = 0; i < numdatastreams; i++)
    setupinfo.antSEFDs[i] = config->getSEFD(i);
  const vector<float> &linesignal = config->getLineSignal();
  for(size_t i = 0; i < LINESIGLEN; i++)
    setupinfo.linesignal[i] = (i < linesignal.size()) ? linesignal[i] : 0;
}

int main(int argc, char* argv[])
{
  // initialize MPI
  MPI_Init(&argc, &argv);

  int numprocs, myid;
  double start, end, elapse;

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  Setup setupinfo;

  // master parse command line argument
  // and send struct setupinfo to all worker processes
  if(myid == MASTER)
  {
    setupinfo.verbose = 0;
    setupinfo.test = 0;
    setupinfo.seed = SEED;
    setupinfo.sfluxdensity = 100;     // source flux density in Jansky
    for(size_t idx = 0; idx < MAXANT; idx++)
      setupinfo.antSEFDs[idx] = 1000;
    for(size_t i = 0; i < LINESIGLEN; i++)
      setupinfo.linesignal[i] = 0;
    setupinfo.numdivs = 1;
    setupinfo.pcal = 0;
    setupinfo.specres = 1;

    // parse command line argument (verbose/test/JSON file name;
    // all other settings come from the JSON configuration)
    cmdparser(argc, argv, setupinfo);
  }

  // MPI_Type_create_struct
  // Create struct for setupinfo
  int block[NITEMS] = {1, 1, 1, 1, MAXANT, MAXLEN, LINESIGLEN, 1, 1, 1};
  MPI_Datatype type[NITEMS] = {MPI_INT, MPI_INT, MPI_UNSIGNED, MPI_FLOAT,
    MPI_INT, MPI_CHAR, MPI_FLOAT, MPI_INT, MPI_INT, MPI_INT};
  MPI_Aint disp[NITEMS];
  MPI_Datatype structtype;

  disp[0] = offsetof(Setup, verbose);
  disp[1] = offsetof(Setup, test);
  disp[2] = offsetof(Setup, seed);
  disp[3] = offsetof(Setup, sfluxdensity);
  disp[4] = offsetof(Setup, antSEFDs);
  disp[5] = offsetof(Setup, inputfilename);
  disp[6] = offsetof(Setup, linesignal);
  disp[7] = offsetof(Setup, numdivs);
  disp[8] = offsetof(Setup, pcal);
  disp[9] = offsetof(Setup, specres);

  MPI_Type_create_struct(NITEMS, block, disp, type, &structtype);
  MPI_Type_commit(&structtype);

  // broadcast setupinfo
  MPI_Bcast(&setupinfo, 1, structtype, MASTER, MPI_COMM_WORLD);

  float tdur = 0.5 * 1e6;
  float testdur = 1; // time duration of test mode
  JsonConfig* config;
  DelayModel* model;
  float dur;                            // duration of the simulated data in seconds
  float durus;                          // duration of the simulated data in microseconds
  float specRes, minStartFreq;
  int numdatastreams, numrecordedbands;
  int numSamps;                         // number of samples to be generated per time block for the common signal
  size_t stime;                         // step time in microsecond == time block
  int configindex = 0;

  double timer = 0.0;
  double tt = 0.0;
  double minprocptrtime = 0.0;
  double refvptime;                     // vptime is the time duration of the samples within the packet
  size_t framespersec;

  // parse the JSON configuration and load the SFXC delay tables
  try
  {
    config = new JsonConfig(setupinfo.inputfilename);
    applyconfig(config, setupinfo);

    vector<string> delaytables;
    for(int i = 0; i < config->getNumDataStreams(); i++)
      delaytables.push_back(config->getDelayTable(i));
    model = new DelayModel(delaytables, config->getStartMJD(),
                           config->getStartSeconds(), config->getScan());
  }
  catch(const exception &e)
  {
    cerr << "Process " << myid << ": " << e.what() << endl;
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    return EXIT_FAILURE;
  }

  // initialize random number generator
  gsl_rng **rng_inst;
  gsl_rng_env_setup();

  rng_inst = (gsl_rng **) malloc(numprocs * sizeof(gsl_rng *));
  rng_inst[myid] = gsl_rng_alloc (gsl_rng_ranlux389);
  gsl_rng_set(rng_inst[myid],setupinfo.seed+myid);

  // retrieve general information from config
  // if in test mode, only generate data for 1 second
  dur = setupinfo.test ? testdur : config->getExecuteSeconds();
  // dur is in second
  // durus is in microsecond
  durus = dur * 1e6;

  // It is possible to use time-based parallelisation
  // Use color to divide the processes into sub-communication groups, where
  // each group generate dur/div seconds of signals for all subbands
  // 'div' is the number of pieces the simulation time is divided into
  int div = setupinfo.numdivs;
  if(myid == MASTER)
  {
    if(numprocs%div != 0)
    {
      cout << "Number of time-based divisions cannot be divided by number of processes ..." << endl;
      MPI_Abort(MPI_COMM_WORLD, ERROR);
    }
    if(div > 1)
    {
      cout << "Use time-based parallelization." << endl;
      cout << "Divide simulation into " << div << " parts" << endl;
    }
  }
  if((size_t)durus%div != 0)
  {
    cout << "WARNING!!!!!!!\n"
         << "Cannot divide durus by div!\n"
         << "Will only generate " << durus/div * div << " microseconds data ..." << endl;
    durus = (float)((size_t)durus/div * div);
  }

  // calculate specRes, number of samples per time block, step time
  if(getSpecRes(config, configindex, specRes, setupinfo.verbose) != EXIT_SUCCESS)
  {
    cout << "Process " << myid << ": Failed to calculate spectral resolution ..." << endl;
    return EXIT_FAILURE;
  };
  // devide the spectral resolution by extra user-defined scaling factor
  specRes /= (float)setupinfo.specres;
  numSamps = getNumSamps(config, configindex, specRes, setupinfo.verbose);
  stime = static_cast<size_t>(1 / specRes); // step time in microsecond
  if(setupinfo.verbose >= 1 && myid == MASTER)
  {
    cout << "SpecRes is " << specRes << " MHz" << endl;
    cout << "number of samples per time block is " << numSamps << "\n"
         << "step time is " << stime << " us" << endl;
    cout << "tdur is " << tdur << endl;
  }
  minStartFreq = getMinStartFreq(config, configindex, setupinfo.verbose);

  // check whether linesignal frequency is out of range
  if(myid == MASTER && setupinfo.linesignal[0] >= numSamps*specRes)
  {
    cout << "ERROR: Line signal frequency value out of range ..." << endl;
    return EXIT_FAILURE;
  };

  numdatastreams = config->getNumDataStreams();

  // use framespersec of the first antenna to calculate vptime
  // use this vptime as time reference
  framespersec = config->getFramesPerSecond(configindex, 0);
  refvptime = 1.0 * 1e6 / framespersec;

  // print out information of the simulation
  // and retrieve information of the number of antenna and subbands
  int sbcount = 0;
  if(myid == MASTER)
  {
    start = MPI_Wtime();
    // array with num-antenna number of elements
    // the value of which is the number of subbands the corresponding antenna has
    vector<int> antsb;
    if(!is_integer(refvptime))
    {
     cout << "VDIF packet time in microsecond is not an integer!! Something is wrong here ..." << endl;
     return EXIT_FAILURE;
    }

    cout << "Generate " << dur << " seconds data for " << numdatastreams << " stations ..." << "\n"
     << "source flux density is " << setupinfo.sfluxdensity << ".\n"
     << "random number generator seed is " << setupinfo.seed << endl;

    for(int i = 0; i < numdatastreams; i++)
    {
     framespersec = (size_t)config->getFramesPerSecond(configindex, i);

     numrecordedbands = config->getDNumRecordedBands(configindex, i);

     cout << "Telescope " << config->getTelescopeName(i) << "\n"
          << " Number of recorded band(s) is " << numrecordedbands << "\n"
          << " Antenna SEFD is " << setupinfo.antSEFDs[i] << "\n"
          << " Number of frames per second is " << framespersec << endl;

     // antenna subbands counter
     int antsbcnt = 0;
     for(int j = 0; j < numrecordedbands; j++)
     {
       sbcount++;
       antsbcnt++;
       cout << "Subband " << j << ":" << "\n"
            << "  Frequency " << config->getDRecordedFreq(configindex, i, j) << "\n"
            << "  Bandwidth " << config->getDRecordedBandwidth(configindex, i, j) << endl;
     }
     antsb.push_back(antsbcnt);
    }
    cout << "Total number of subbands of all antennas is " << sbcount << endl;

    // general information from model
    if(setupinfo.verbose >= 2)
    {
     cout << "Total number of subbands is  " << sbcount << endl;
     cout << "Number of scans is " << model->getNumScans() << endl;
     cout << "Scan start second is " << model->getScanStartSec(0, config->getStartMJD(), config->getStartSeconds()) << endl;
     cout << "Scan end second is " << model->getScanEndSec(0, config->getStartMJD(), config->getStartSeconds()) << endl;
    }
  }
  // Broadcast sbcount to all processes
  MPI_Bcast(&sbcount, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

  durus /= (float)div;
  int color = myid / (numprocs/div);
  // Create communication groups for each time-based partition
  MPI_Comm local_comm;
  MPI_Comm_split(MPI_COMM_WORLD, color, myid, &local_comm);
  int local_rank, local_size;
  MPI_Comm_rank(local_comm, &local_rank);
  MPI_Comm_size(local_comm, &local_size);
  if(myid == MASTER && sbcount%local_size != 0)
  {
    cout << "ERROR: number of processes cannot be divided by total nubmer of subbands ...\n"
         << "If time-based serialisation is specified, please make sure np/numdivs can be divided by total number of subbands"<< endl;
    return EXIT_FAILURE;
  }

  int numsbperproc = sbcount / local_size;
  if(setupinfo.verbose >= 1)
    cout << "number of subbands per process is " << numsbperproc << endl;
  // subbandsinfo is only used by MASTER to store subband information
  // and distribute the information to all processes using MPI_Scatter(...)

  // create a 2D array with the totoal number of subbands
  // containting antenna index and subband index of each subband
  int* subbandsinfo;
  int* sbinfo = new int [numsbperproc * 2];
  if(local_rank == MASTER)
  {
    subbandsinfo = new int [sbcount * 2];

    int numprocessed = 0;
    for(int i = 0; i < numdatastreams; i++)
    {
      numrecordedbands = config->getDNumRecordedBands(configindex, i);
      for(int j = 0; j < numrecordedbands; j ++)
      {
        int idx = 2 * (numprocessed + j);
        subbandsinfo[idx] = i;
        subbandsinfo[idx+1] = j;
      }
      numprocessed += numrecordedbands;
    }
  }

  // Scatter sbinfo to each process
  // Disbribute array of (antidx, sbidx) information to each process
  MPI_Scatter(&subbandsinfo[0], numsbperproc*2, MPI_INT, &sbinfo[0], numsbperproc*2, MPI_INT, MASTER, local_comm);

  // master generates common signal
  // each process copies its corresponding subband data to the right place
  // if use time-based parallelization
  // master is local_rank 0 of each group

  size_t stdur = tdur/stime;

  // Subband-based parallelization starts here
  vector<Subband*> subbands;

  // every process initialize its own subbands
  initSubbands(config, configindex, model, specRes, minStartFreq, subbands, numsbperproc, tdur, setupinfo, sbinfo, color, durus);

  if(setupinfo.verbose >= 1)
    cout << "Process " << myid
         << ", local rank " << local_rank << " of group " << color << endl;

  float* commFreqSig;                  // 0.5 seconds common frequency domain signal
  float* commSlice;

  // allocate memory for the common frequency domain signal
  int sampsize = numSamps*2*stdur;          // *2 due to complex number
  commFreqSig = new float [sampsize];
  commSlice = new float [numSamps*2];

  // apply Gaussian filter as spectral line
  float* gaussianfilter = new float [numSamps*2];
  if(setupinfo.linesignal[0] > EPSILON)
  {
    gengaussianfilter(gaussianfilter, setupinfo.linesignal, numSamps, specRes);
  }

  while(timer < durus)
  {
    if(local_rank == MASTER)
    {
      if(setupinfo.verbose >= 1)
        cout << "Generate " << tdur << " us signal" << endl;
      gencplx(commFreqSig, sampsize, STDEV, rng_inst[myid], setupinfo.verbose);

      if(setupinfo.linesignal[0] > EPSILON)
      {
        for(size_t t = 0; t < (size_t)sampsize/numSamps/2; t++)
          for(size_t idx = 0; idx < (size_t)numSamps*2; idx++)
            commFreqSig[t*numSamps*2+idx] *= gaussianfilter[idx];
      }

    }

    MPI_Bcast(commFreqSig, sampsize, MPI_FLOAT, MASTER, MPI_COMM_WORLD);

    vector<Subband*>::iterator it;
    for(it = subbands.begin(); it != subbands.end(); ++it)
    {
      for(size_t t = 0; t < stdur; t++)
      {
        for(size_t samp = 0; samp < (size_t)numSamps*2; samp++)
        {
          size_t idx = t*numSamps*2+samp;
          commSlice[samp] = commFreqSig[idx];
        }
        (*it)->fabricatedata(commSlice, rng_inst[myid], setupinfo.sfluxdensity);
      }
      // after TDUR time signal is generated for each subband array
      // set the current pointer of each array back to the beginning of the second half

      (*it)->setcptr((*it)->getlength() / 2);
      if(setupinfo.verbose >= 2)
        cout << "Antenna " << (*it)->getantIdx() << " subband " << (*it)->getsbIdx()
                           << " set current pointer back to " << (*it)->getlength() / 2 << endl;
    }
    // move data in each array from the second half to the first half
    // and set the process pointer to the proper location
    // i.e. data is moved half array ahead, therefore process pointer
    // is moved half array ahead
    movedata(subbands, setupinfo.verbose);

    // each subband calculates its own process pointer time
    minprocptrtime = getMinProcPtrTime(subbands, setupinfo.verbose);
    if(setupinfo.verbose >= 2) cout << "Process " << local_rank << ": local minimum procptrtime is " << minprocptrtime << endl;

    MPI_Allreduce(&minprocptrtime, &tt, 1, MPI_DOUBLE, MPI_MIN, local_comm);
    if(setupinfo.verbose >= 2)
      cout << "Process " << local_rank <<": tt is " << tt << ", timer is : " << timer
           << ", tdur is " << tdur << ", durus is " << durus << endl;

    if((local_rank == MASTER) && (tt >= tdur))
    {
      cout << "**the lowest process pointer time is larger than tdur!!!\n"
              "**Something is wrong here!!!" << endl;
      MPI_Abort(MPI_COMM_WORLD, ERROR);
      return (EXIT_FAILURE);
    }


    while((tt < tdur) && (timer < durus))
    {
      // process and packetize one vdif packet for each subband array
      int rc = processAndPacketize(subbands, model, setupinfo.verbose, setupinfo.pcal);
      if(rc)
      {
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        return(EXIT_FAILURE);
      }
      tt += refvptime;
      timer += refvptime;
      if((local_rank == MASTER) && (setupinfo.verbose >=2))
        cout << "tt is " << tt << ", timer is " << timer << endl;
    }

    MPI_Bcast(&timer, 1, MPI_DOUBLE, MASTER, local_comm);
  }


  // free allocated common signal memory
  delete commFreqSig;

  // combine VDIF files
  if(local_size < numdatastreams)
  {
    size_t antperproc = (numdatastreams%local_size == 0) ? (size_t)numdatastreams/local_size : (size_t)numdatastreams/local_size + 1;
    for(size_t idx = 0; idx < antperproc; idx++)
    {
      size_t antidx = idx + local_rank * antperproc;
      if(antidx < (size_t)numdatastreams)
        vdifzipper(config, configindex, durus, setupinfo.verbose, antidx, color);
    }
  }
  else
  {
    if(local_rank < numdatastreams)
      vdifzipper(config, configindex, durus, setupinfo.verbose, local_rank, color);
  }

  //if(local_rank < numdatastreams)
  //  vdifzipper(config, configindex, durus, setupinfo.verbose, local_rank, color);

  MPI_Type_free(&structtype);
  MPI_Comm_free(&local_comm);

  if(numprocs < numdatastreams)
  {
    size_t antperproc = (numdatastreams%numprocs == 0) ? (size_t)numdatastreams/numprocs : (size_t)numdatastreams/numprocs + 1;
    if(setupinfo.verbose >= 1) cout << "catvdif: number of antennas to process for process " << myid << " is "<< antperproc << endl;
    for(size_t idx = 0; idx < antperproc; idx++)
    {
      size_t antidx = idx + myid * antperproc;
      if(antidx < (size_t)numdatastreams)
      {
        string antname = config->getTelescopeName(antidx);
        catvdif(antname, setupinfo.verbose, antidx, div);
      }
    }
  }
  else
  {
    if(myid < numdatastreams)
    {
      string antname = config->getTelescopeName(myid);
      catvdif(antname, setupinfo.verbose, myid, div);
    }
  }

  //if(myid < numdatastreams)
  //  catvdif(config, configindex, durus, setupinfo.verbose, myid, div);

  freeSubbands(subbands);
  delete model;
  delete config;
  if(myid == MASTER)
  {
    end = MPI_Wtime();
    elapse = (end - start)/60.0; // convert time from seconds to minutes

    cout << "All data has been generated successfully, bye!" << endl;
    cout << "Total duration is " << elapse << " minutes." << endl;
  }

  // free random number generator
  gsl_rng_free(rng_inst[myid]);
  MPI_Finalize();

  return(EXIT_SUCCESS);
}
