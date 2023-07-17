#include "HitSelectorSpace.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <IMPL/LCCollectionVec.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTrackerConf.h>

#include "TMath.h"
#include "TVector3.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

HitSelectorSpace aHitSelectorSpace;

HitSelectorSpace::HitSelectorSpace() : Processor("HitSelectorSpace")
{

    // Modify processor description
    _description = "HitSelectorSpace applies space selections to reduce the BIB";

    // Input collection
    registerProcessorParameter("TrackerHitCollectionName",
                               "Name of the TrackerHit input collection",
                               m_inputHitCollection,
                               std::string("VertexBarrelCollection"));

    // Output collection
    registerProcessorParameter("GoodHitCollection",
                               "Good hits from tracker",
                               m_outputHitCollection,
                               std::string("VertexBarrelGoodCollection"));
}

void HitSelectorSpace::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;
}

void HitSelectorSpace::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void HitSelectorSpace::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG8) << "Processing event " << _nEvt << std::endl;

    // Get the collection of tracker hits
    LCCollection *trackerHitCollection = 0;
    getCollection(trackerHitCollection, m_inputHitCollection, evt);

    std::string encoderString = trackerHitCollection->getParameters().getStringVal("CellIDEncoding");
    UTIL::CellIDDecoder<TrackerHitPlane> myCellIDEncoding(encoderString);

    // Make the output collections
    LCCollectionVec *GoodHitsCollection = new LCCollectionVec(trackerHitCollection->getTypeName());
    GoodHitsCollection->setSubset(true);
    GoodHitsCollection->parameters().setValue("CellIDEncoding", encoderString);

    int nHits = trackerHitCollection->getNumberOfElements();
    // Set the map of responses 
    memset(&m_accepted, false, nHits);

    // First sort hits in a map
    m_hitsMap.clear();
    for (int itHit = 0; itHit < nHits; itHit++)
    {
        // Get the hit
        TrackerHitPlane *hit = static_cast<TrackerHitPlane *>(trackerHitCollection->getElementAt(itHit));
        
        unsigned int layer= myCellIDEncoding(hit)["layer"];
        unsigned int side= myCellIDEncoding(hit)["side"];
        unsigned int ladder= myCellIDEncoding(hit)["module"];
        unsigned int module= myCellIDEncoding(hit)["sensor"];

        MySensorPos sensPos = {layer, side, ladder, module};
        if (m_hitsMap.find(sensPos) == m_hitsMap.end())
        {
            m_hitsMap[sensPos] = std::vector<size_t>();
            m_hitsMap[sensPos].reserve(nHits);
        }
        m_hitsMap[sensPos].push_back(itHit);

    }

    // Loop over tracker hits
    for (int itHit = 0; itHit < nHits; itHit++)
    {

        // Skip accepted hits
        if (m_accepted[itHit]){
            streamlog_out(DEBUG0) << "Skipping already accepted hit" << std::endl;
            continue;
        }

        // Get the hit
        TrackerHitPlane *hit = static_cast<TrackerHitPlane *>(trackerHitCollection->getElementAt(itHit));
        unsigned int layer = myCellIDEncoding(hit)["layer"];
        unsigned int side = myCellIDEncoding(hit)["side"];
        unsigned int ladder = myCellIDEncoding(hit)["module"];
        unsigned int module = myCellIDEncoding(hit)["sensor"];

        streamlog_out(DEBUG0) << "Hit position " << layer << " " << ladder << " " << module << std::endl;

        // We go inside out and skip the outer layers
        if(layer==1 || layer==3 || layer==5 ||layer==7){
            streamlog_out(DEBUG0) << "Skipping hit in outer layer of pair" << std::endl;
            continue;
        }

        // get the hit position
        TVector3 pos(hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);
        double min_dR = 999999.;
        double dtheta_closest = 0.;
        double dphi_closest = 0.;
        double theta = pos.Theta();
        double dtheta_cut = 0.01;
        double dphi_cut = 0.001;

        unsigned int the_other_layer = 1;
        if (layer == 2)
        {
            the_other_layer = 3;
            dtheta_cut = 0.005;
            dphi_cut = 0.001;
        }
        if (layer == 4)
        {
            the_other_layer = 5;
            dtheta_cut = 0.002;
            dphi_cut = 0.001;
        }
        if (layer == 6)
        {
            the_other_layer = 7;
            dtheta_cut = 0.001;
            dphi_cut = 0.001;
        }

        const MySensorPos theOtherPos = {the_other_layer, side, ladder, module};

        // Checking if there are any hits in the other layer
        if (m_hitsMap.find(theOtherPos) == m_hitsMap.end())
        {
            streamlog_out(DEBUG0) << "No hits in outer layer of pair" << std::endl;
            continue;
        }

        for (size_t jitHit : m_hitsMap.at(theOtherPos))
        {
            TrackerHitPlane *hit2 = static_cast<TrackerHitPlane *>(trackerHitCollection->getElementAt(jitHit));
            TVector3 pos2(hit2->getPosition()[0], hit2->getPosition()[1], hit2->getPosition()[2]);
            double dtheta = pos2.Theta() - theta;
            double dphi = pos.DeltaPhi(pos2);
            double dR = sqrt(dphi * dphi + dtheta * dtheta);
            if (dR < min_dR)
            {
                min_dR = dR;
                dtheta_closest = dtheta;
                dphi_closest = dphi;
            }
            if (fabs(dtheta) > dtheta_cut){
                streamlog_out(DEBUG0) << " -> fail dtheta " << fabs(dtheta) << std::endl;
                continue;
            }
            if (fabs(dphi) > dphi_cut){
                streamlog_out(DEBUG0) << " -> fail dphi " << fabs(dphi) << std::endl;
                continue;
            }
            unsigned int layer2 = myCellIDEncoding(hit2)["layer"];

            streamlog_out(DEBUG0) << " --> accepted hit in outer layer (" << layer2 << ") of pair with " << dtheta << " " << dphi << std::endl;
            m_accepted[jitHit] = true;
        }

        if (fabs(dtheta_closest) < dtheta_cut && fabs(dphi_closest) < dphi_cut)
        {
            streamlog_out(DEBUG0) << " --> accepted hit in inner layer of pair with " << fabs(dtheta_closest) << std::endl;
            m_accepted[itHit] = true;
        }

    }

    // Once more to add the hits to the output (in slices)
    for (size_t itHit = 0; itHit < nHits; itHit++)
    {
        TrackerHitPlane *hit = static_cast<TrackerHitPlane *>(trackerHitCollection->getElementAt(itHit));

        if (m_accepted[itHit])
        {
            GoodHitsCollection->addElement(hit);
        }
    }

    // Store the filtered hit collections
    evt->addCollection(GoodHitsCollection, m_outputHitCollection);

    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    streamlog_out(DEBUG) << "   done processing event: " << evt->getEventNumber()
                         << "   in run:  " << evt->getRunNumber() << std::endl;
    
    _nEvt++;
}

void HitSelectorSpace::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void HitSelectorSpace::end()
{

    //   std::cout << "HitSelectorSpace::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void HitSelectorSpace::getCollection(LCCollection *&collection, std::string collectionName, LCEvent *evt)
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
