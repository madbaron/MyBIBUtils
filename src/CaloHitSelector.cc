#include "CaloHitSelector.h"
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <filesystem>

#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>

#include <UTIL/LCRelationNavigator.h>
#include <IMPL/LCCollectionVec.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTrackerConf.h>
#include <IMPL/LCRelationImpl.h>

#include "TH1D.h"
#include "TVector3.h"
#include "TMath.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

CaloHitSelector aCaloHitSelector;

CaloHitSelector::CaloHitSelector() : Processor("CaloHitSelector")
{

    // Modify processor description
    _description = "CaloHitSelector applies E selections to reduce the BIB";

    // Input collection
    registerProcessorParameter("CaloHitCollectionName",
                               "Name of the CalorimeterHit input collection",
                               m_inputHitCollection,
                               std::string("EcalBarrelCollectionRec"));

    // Output collection
    registerProcessorParameter("GoodHitCollection",
                               "Good hits from calo",
                               m_outputHitCollection,
                               std::string("EcalBarrelCollectionSel"));

    // Input relation collection
    registerProcessorParameter("CaloRelationCollectionName",
                               "Name of the CalorimeterHit input relation collection",
                               m_inputRelationCollection,
                               std::string("EcalBarrelRelationsSimRec"));

    // Output relation collection
    registerProcessorParameter("GoodRelationCollection",
                               "Good hits SimRec relations",
                               m_outputRelationCollection,
                               std::string("EcalBarrelRelationsSimSel"));

    // ROOT map
    registerProcessorParameter("ThresholdsFilePath",
                               "Path to ROOT file",
                               m_thFile,
                               std::string(""));

    // N sigma for dynamic threshold
    registerProcessorParameter("Nsigma",
                               "Number of BIB E sigma",
                               m_Nsigma,
                               3);

    // Fixed threshold
    registerProcessorParameter("FlatThreshold",
                               "Cut in GeV",
                               m_FlatThreshold,
                               0.);

    //TimeWindowMin
    registerProcessorParameter("TimeWindowMin",
                               "Minimum time window for hit selection",
                               m_time_windowMin,
                               -0.5);

    //TimeWindowMax
    registerProcessorParameter("TimeWindowMax",
                               "Maximum time window for hit selection",
                               m_time_windowMax,
                               10.);

    // Subtract expected BIB energy
    registerProcessorParameter("DoBIBsubtraction",
                               "Correct cell energy for mean expected BIB contribution",
                               m_doBIBsubtraction,
                               bool(false));
}

void CaloHitSelector::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;

    // open ROOT file and get threshold histograms
    m_th_file = new TFile(m_thFile.c_str());
    m_thresholdMap = (TH2D *)m_th_file->Get("th_2dmode_sym");
    m_thresholdMap->SetDirectory(0);
    m_stddevMap = (TH2D *)m_th_file->Get("stddev_sym");
    m_stddevMap->SetDirectory(0);
    m_th_file->Close();
}

void CaloHitSelector::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void CaloHitSelector::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG) << "Processing event " << _nEvt << std::endl;
    streamlog_out(DEBUG) << " in " << this->name() << std::endl;

    // Get the collection of calo hits
    LCCollection *caloHitCollection = 0;
    getCollection(caloHitCollection, m_inputHitCollection, evt);

    LCCollection *inputHitRel = 0;
    getCollection(inputHitRel, m_inputRelationCollection, evt);

    LCCollectionVec *outputHitCol = 0;
    LCCollection *outputHitRel = 0;

    if (caloHitCollection != 0 && inputHitRel != 0)
    {

        std::string encoderString = caloHitCollection->getParameters().getStringVal(LCIO::CellIDEncoding);
        UTIL::CellIDDecoder<CalorimeterHit> myCellIDEncoding(encoderString);

        // Make the output collections
        outputHitCol = new LCCollectionVec(caloHitCollection->getTypeName());
        outputHitCol->setSubset(true);
        outputHitCol->parameters().setValue(LCIO::CellIDEncoding, encoderString);
        outputHitCol->setFlag(outputHitCol->getFlag() | (1 << EVENT::LCIO::CHBIT_LONG));
        outputHitCol->setFlag(outputHitCol->getFlag() | (1 << EVENT::LCIO::RCHBIT_TIME));

        // reco-sim relation output collections
        UTIL::LCRelationNavigator thitNav = UTIL::LCRelationNavigator( LCIO::CALORIMETERHIT, LCIO::SIMCALORIMETERHIT );

        int nHits = caloHitCollection->getNumberOfElements();

        // Now loop over hits again applying threshold
        for (int itHit = 0; itHit < nHits; itHit++)
        {
            // Get the hit
            CalorimeterHit *hit = static_cast<CalorimeterHit *>(caloHitCollection->getElementAt(itHit));
            unsigned int layer = myCellIDEncoding(hit)["layer"];

            // hit position
            TVector3 hitPos(hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);
            double hit_theta = hitPos.Theta();
            if (hit_theta > TMath::Pi() / 2) // map is symmetrized around pi/2
            {
                hit_theta = TMath::Pi() - hit_theta;
            }

            unsigned int binx = m_thresholdMap->GetXaxis()->FindBin(hit_theta);
            unsigned int biny = m_thresholdMap->GetYaxis()->FindBin(layer);

            double threshold = m_thresholdMap->GetBinContent(binx, biny) + m_Nsigma * m_stddevMap->GetBinContent(binx, biny);
            if (m_FlatThreshold > 0.)
            {
                threshold = m_FlatThreshold;
            }

            double correction = m_thresholdMap->GetBinContent(binx, biny);

            double hit_energy = hit->getEnergy();
            if (m_doBIBsubtraction)
            {
                hit_energy = hit_energy - correction;
            }

            if (hit_energy > threshold)
            {
                // Compute time correction
                float timeCorrection(0);
                float r(0);
                for (int i=0; i<3; i++)
                    r+=pow(hit->getPosition()[i],2);
                timeCorrection = sqrt(r)/TMath::C(); // [speed of light in mm/ns]

                float relativetime = hit->getTime() - timeCorrection; // wrt time of flight

                if (relativetime>m_time_windowMin && relativetime<m_time_windowMax){

                    streamlog_out(DEBUG0) << " accepted hit " << hit_energy << " theta " << hit_theta << std::endl;

                    outputHitCol->addElement(hit);

                    LCRelation *rel = static_cast<LCRelation *>(inputHitRel->getElementAt(itHit));
                    SimCalorimeterHit *simhit = static_cast<SimCalorimeterHit *>(rel->getTo());

                    thitNav.addRelation(hit, simhit);
                }
            }
        }

        // Store the filtered hit collections
        evt->addCollection(outputHitCol, m_outputHitCollection);
        outputHitRel = thitNav.createLCCollection();
        evt->addCollection(outputHitRel, m_outputRelationCollection);
    }

    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    streamlog_out(DEBUG) << "   done processing event: " << evt->getEventNumber()
                         << "   in run:  " << evt->getRunNumber() << std::endl;

    _nEvt++;
}

void CaloHitSelector::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void CaloHitSelector::end()
{
    //   std::cout << "CaloHitSelector::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void CaloHitSelector::getCollection(LCCollection *&collection, const std::string &collectionName, LCEvent *evt)
{
    try
    {
        collection = evt->getCollection(collectionName);
    }
    catch (DataNotAvailableException &e)
    {
        streamlog_out(DEBUG) << "- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
        return;
    }
    return;
}
