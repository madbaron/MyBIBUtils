#include "CaloHitSelector.h"
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <filesystem>

#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/CalorimeterHitImpl.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTrackerConf.h>
#include <IMPL/LCRelationImpl.h>

#include "TH1D.h"
#include "TVector3.h"

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

    if (caloHitCollection != 0 && inputHitRel != 0)
    {

        std::string encoderString = caloHitCollection->getParameters().getStringVal("CellIDEncoding");
        UTIL::CellIDDecoder<CalorimeterHit> myCellIDEncoding(encoderString);

        // Make the output collections
        LCCollectionVec *GoodHitsCollection = new LCCollectionVec("CalorimeterHit");
        GoodHitsCollection->setSubset(false);
        GoodHitsCollection->parameters().setValue("CellIDEncoding", encoderString);

        // reco-sim relation output collections
        LCCollectionVec *outputHitRel = new LCCollectionVec("LCRelation");
        outputHitRel->parameters().setValue("FromType", "CalorimeterHit");
        outputHitRel->parameters().setValue("ToType", "SimCalorimeterHit");
        LCFlagImpl lcFlag_rel(inputHitRel->getFlag());
        outputHitRel->setFlag(lcFlag_rel.getFlag());

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
            double correction = m_thresholdMap->GetBinContent(binx, biny);

            if (hit->getEnergy() > threshold)
            {
                streamlog_out(DEBUG0) << " accepted hit " << hit->getEnergy() << " theta " << hit_theta << std::endl;

                CalorimeterHitImpl *hit_new = new CalorimeterHitImpl();

                hit_new->setCellID0(hit->getCellID0());
                hit_new->setCellID1(hit->getCellID1());
                hit_new->setType(hit->getType());
                hit_new->setRawHit(hit->getRawHit());
                hit_new->setPosition(hit->getPosition());
                hit_new->setTime(hit->getTime());
                hit_new->setEnergy(hit->getEnergy() - correction);
                hit_new->setEnergyError(hit->getEnergyError());

                GoodHitsCollection->addElement(hit_new);

                LCRelation *rel = dynamic_cast<LCRelation *>(inputHitRel->getElementAt(itHit));
                LCRelationImpl *rel_new = new LCRelationImpl();

                rel_new->setFrom(hit_new);
                rel_new->setTo(rel->getTo());
                rel_new->setWeight(rel->getWeight());

                outputHitRel->addElement(rel_new);
            }
        }

        // Store the filtered hit collections
        evt->addCollection(GoodHitsCollection, m_outputHitCollection);
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
