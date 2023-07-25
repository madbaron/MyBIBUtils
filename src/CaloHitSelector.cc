#include "CaloHitSelector.h"
#include <iostream>
#include <vector>
#include <map>
#include <math.h>

#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/CalorimeterHitImpl.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTrackerConf.h>
#include <IMPL/LCRelationImpl.h>

#include <marlin/AIDAProcessor.h>

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

    // N calo layers
    registerProcessorParameter("Nlayers",
                               "Number of calorimeter layers",
                               m_Nlayer,
                               50);

    // N sigma for dynamic threshold
    registerProcessorParameter("Nsigma",
                               "Number of BIB E sigma",
                               m_Nsigma,
                               3);

    // Save histograms
    registerProcessorParameter("SaveHistograms",
                               "Flag to save the threshold histograms",
                               m_fillHistos,
                               false);
}

void CaloHitSelector::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;

    // --- Initialize the AIDAProcessor and book the diagnostic histograms:

    AIDAProcessor::histogramFactory(this);

    size_t nThetaBins = arrBins_theta.size() - 1;
    m_thresholdMap = new TH2D("Threshold_vs_theta_layer", "Threshold_vs_theta_layer", nThetaBins, 0, nThetaBins, m_Nlayer, 0, m_Nlayer);
    m_correctionMap = new TH2D("Correction_vs_theta_layer", "Correction_vs_theta_layer", nThetaBins, 0, nThetaBins, m_Nlayer, 0, m_Nlayer);
}

void CaloHitSelector::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void CaloHitSelector::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG) << "Processing event " << _nEvt << std::endl;

    // Get the collection of calo hits
    LCCollection *caloHitCollection = 0;
    getCollection(caloHitCollection, m_inputHitCollection, evt);

    LCCollection *inputHitRel = 0;
    getCollection(inputHitRel, m_inputRelationCollection, evt);

    std::string encoderString = caloHitCollection->getParameters().getStringVal("CellIDEncoding");
    UTIL::CellIDDecoder<CalorimeterHit> myCellIDEncoding(encoderString);

    // Make the output collections
    LCCollectionVec *GoodHitsCollection = new LCCollectionVec(caloHitCollection->getTypeName());
    GoodHitsCollection->setSubset(false);
    GoodHitsCollection->parameters().setValue("CellIDEncoding", encoderString);

    // reco-sim relation output collections
    LCCollectionVec *outputHitRel = new LCCollectionVec(inputHitRel->getTypeName());
    outputHitRel->parameters().setValue("FromType", inputHitRel->parameters().getStringVal("FromType"));
    outputHitRel->parameters().setValue("ToType", inputHitRel->parameters().getStringVal("ToType"));
    LCFlagImpl lcFlag_rel(inputHitRel->getFlag());
    outputHitRel->setFlag(lcFlag_rel.getFlag());

    std::vector<TH2D *> layer_map;
    for (int iLayer = 0; iLayer < m_Nlayer; iLayer++)
    {
        TString name;
        name.Form("%s%s%i", m_inputHitCollection, "hit_E_vs_theta_layer_", iLayer);
        TH2D *histo = new TH2D(name, name, arrBins_theta.size() - 1, arrBins_theta.data(), 100, 0, 0.1);
        layer_map.push_back(histo);
    }

    // Loop over calo hits
    int nHits = caloHitCollection->getNumberOfElements();
    for (int itHit = 0; itHit < nHits; itHit++)
    {
        // Get the hit
        CalorimeterHit *hit = static_cast<CalorimeterHit *>(caloHitCollection->getElementAt(itHit));
        unsigned int layer = myCellIDEncoding(hit)["layer"];

        streamlog_out(DEBUG0) << " " << std::endl;
        streamlog_out(DEBUG0) << " Found hit L " << layer << std::endl;

        // hit position
        TVector3 hitPos(hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);
        double hit_theta = hitPos.Theta();

        streamlog_out(DEBUG0) << " E " << hit->getEnergy() << " theta " << hit_theta << std::endl;

        layer_map[layer]->Fill(hit_theta, hit->getEnergy());
    }

    // Now extract thresholds
    std::vector<std::vector<double>> threshold_map;
    std::vector<std::vector<double>> correction_map;
    for (int iLayer = 0; iLayer < m_Nlayer; iLayer++)
    {
        std::vector<double> threshold_map_theta;
        std::vector<double> correction_map_theta;
        size_t maxBin = arrBins_theta.size() - 1;
        for (size_t iBin = 0; iBin < maxBin; iBin++)
        {
            TString name_proj;
            name_proj.Form("%s%li", "_py_", iBin);
            TH1D *h_my_proj = layer_map[iLayer]->ProjectionY(name_proj, iBin, iBin + 1);
            threshold_map_theta.push_back(h_my_proj->GetMean() + m_Nsigma * h_my_proj->GetStdDev());
            correction_map_theta.push_back(h_my_proj->GetMean());
            delete h_my_proj;
        }
        threshold_map.push_back(threshold_map_theta);
        correction_map.push_back(correction_map_theta);
    }

    // Now loop over hits again applying threshold
    for (int itHit = 0; itHit < nHits; itHit++)
    {
        // Get the hit
        CalorimeterHit *hit = static_cast<CalorimeterHit *>(caloHitCollection->getElementAt(itHit));
        unsigned int layer = myCellIDEncoding(hit)["layer"];

        // hit position
        TVector3 hitPos(hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);
        double hit_theta = hitPos.Theta();

        unsigned int binx = layer_map[layer]->GetXaxis()->FindBin(hit_theta);

        if (hit->getEnergy() > threshold_map[layer][binx - 1])
        {
            streamlog_out(DEBUG0) << " accepted hit " << hit->getEnergy() << " theta " << hit_theta << std::endl;

            CalorimeterHitImpl *hit_new = new CalorimeterHitImpl();

            hit_new->setCellID0(hit->getCellID0());
            hit_new->setCellID1(hit->getCellID1());
            hit_new->setType(hit->getType());
            hit_new->setRawHit(hit->getRawHit());
            hit_new->setPosition(hit->getPosition());
            hit_new->setTime(hit->getTime());
            hit_new->setEnergy(hit->getEnergy() - correction_map[layer][binx - 1]);
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

    if (m_fillHistos)
    {
        for (int iLayer = 0; iLayer < m_Nlayer; iLayer++)
        {
            size_t maxBin = arrBins_theta.size() - 1;
            for (size_t iBin = 0; iBin < maxBin; iBin++)
            {
                m_thresholdMap->Fill(iBin, iLayer, threshold_map[iLayer][iBin]);
                m_correctionMap->Fill(iBin, iLayer, correction_map[iLayer][iBin]);
            }
        }
    }

    // Store the filtered hit collections
    evt->addCollection(GoodHitsCollection, m_outputHitCollection);
    evt->addCollection(outputHitRel, m_outputRelationCollection);

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
        streamlog_out(DEBUG5) << "- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
        return;
    }
    return;
}
