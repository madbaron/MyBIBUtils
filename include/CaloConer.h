#ifndef CaloConer_h
#define CaloConer_h 1

#include "marlin/Processor.h"
#include "lcio.h"

#include <string>
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"

using namespace lcio;
using namespace marlin;

/**  Example processor for marlin.
 *
 *  If compiled with MARLIN_USE_AIDA
 *  it creates a histogram (cloud) of the MCParticle energies.
 *
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4>
 *  A histogram.
 *
 * @param CollectionName Name of the MCParticle collection
 *
 * @author F. Meloni, DESY;
 * @version $Id: CaloConer.h,v 0.1 2020-09-27 11:24:21 fmeloni Exp $
 */

class CaloConer : public Processor
{

public:
  virtual Processor *newProcessor() { return new CaloConer; }

  CaloConer();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

  // Call to get collections
  void getCollection(LCCollection *&, const std::string &, LCEvent *);

protected:
  // Collection names for (in/out)put
  std::string m_inputMCParticleCollection = "";
  std::string m_inputHitCollection = "";
  std::string m_outputHitCollection = "";
  std::string m_inputRelationCollection = "";
  std::string m_outputRelationCollection = "";

  double m_ConeSize = 0.2;

  int _nRun{};
  int _nEvt{};
};

#endif
