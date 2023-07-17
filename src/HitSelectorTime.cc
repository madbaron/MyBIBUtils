#include "HitSelectorTime.h"
#include <iostream>
#include "TMath.h"

#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>

#include <IMPL/LCCollectionVec.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTrackerConf.h>

#include "TMath.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

HitSelectorTime aHitSelectorTime;

HitSelectorTime::HitSelectorTime() : Processor("HitSelectorTime")
{

    // Modify processor description
    _description = "HitSelectorTime applies time selections to reduce the BIB";

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

void HitSelectorTime::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;
}

void HitSelectorTime::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void HitSelectorTime::processEvent(LCEvent *evt)
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

    // Loop over tracker hits
    int nHits = trackerHitCollection->getNumberOfElements();
    for (int itHit = 0; itHit < nHits; itHit++)
    {
        // Get the hit
        TrackerHitPlane *hit = static_cast<TrackerHitPlane *>(trackerHitCollection->getElementAt(itHit));
        unsigned int layer = myCellIDEncoding(hit)["layer"];
        unsigned int subdet = myCellIDEncoding(hit)["system"];
        unsigned int module = myCellIDEncoding(hit)["module"];
        unsigned int side = myCellIDEncoding(hit)["side"];
        unsigned int sensor = myCellIDEncoding(hit)["sensor"];

        streamlog_out(DEBUG0) << " " << std::endl;
        streamlog_out(DEBUG0) << " Found hit L " << layer << " Su " << subdet << " M " << module << " Si " << side << " Se " << sensor << std::endl;

        // hit position
        double r = sqrt(hit->getPosition()[0] * hit->getPosition()[0] + hit->getPosition()[1] * hit->getPosition()[1]);
        streamlog_out(DEBUG0) << " E " << hit->getEDep() << " time " << hit->getTime() << " r " << r << std::endl;
        double t_fly = r * 1.E6 / TMath::C();
        double t_arr = hit->getTime() - t_fly + 0.2167; // ugly should implement in digitizer

        streamlog_out(DEBUG0) << " t " << hit->getTime() << " t_fly " << t_fly << " t_arr " << t_arr << std::endl;

        // The following conditions are tentative
        if (t_arr > -0.15 && t_arr < 0.15){
            streamlog_out(DEBUG0) << " --> accepted" << std::endl;
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

void HitSelectorTime::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void HitSelectorTime::end()
{

    //   std::cout << "HitSelectorTime::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void HitSelectorTime::getCollection(LCCollection *&collection, std::string collectionName, LCEvent *evt)
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
