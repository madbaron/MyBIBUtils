#include "HitSlimmer.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>

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

HitSlimmer aHitSlimmer;

HitSlimmer::HitSlimmer() : Processor("HitSlimmer")
{

    // Modify processor description
    _description = "HitSlimmer outputs collection of unused hits";

    // Input collection
    registerProcessorParameter("HitsCollectionName",
                               "Name of the Hits input collection",
                               m_inputHitCollection,
                               std::string("HitsCollection"));

    registerProcessorParameter("TrackCollectionName",
                               "Name of reconstructed track input collection",
                               m_inputTrackCollection,
                               std::string("Tracks"));

    // Output collection
    registerProcessorParameter("SlimmedHitsCollectionName",
                               "Name of the slimmed hits output collection",
                               m_outputHitCollection,
                               std::string("SlimmedHits"));
}

void HitSlimmer::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;
}

void HitSlimmer::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void HitSlimmer::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG8) << "Processing event " << _nEvt << std::endl;

    // Get the collection of tracker hits
    LCCollection *trackerHitCollection = 0;
    getCollection(trackerHitCollection, m_inputHitCollection, evt);

    std::string encoderString = trackerHitCollection->getParameters().getStringVal("CellIDEncoding");
    UTIL::CellIDDecoder<TrackerHitPlane> myCellIDEncoding(encoderString);

    // Get the collection of tracks
    LCCollection *trackCollection = 0;
    getCollection(trackCollection, m_inputTrackCollection, evt);

    // vector of used hits
    LCCollectionVec *UsedHitsCollection = new LCCollectionVec(trackerHitCollection->getTypeName());
    UsedHitsCollection->parameters().setValue("CellIDEncoding", encoderString);
    UsedHitsCollection->setTransient(true);

    // Make the output collections
    LCCollectionVec *SlimmedHitsCollection = new LCCollectionVec(trackerHitCollection->getTypeName());
    SlimmedHitsCollection->setSubset(true);
    SlimmedHitsCollection->parameters().setValue("CellIDEncoding", encoderString);

    int nTracks = trackCollection->getNumberOfElements();
    streamlog_out(DEBUG) << "  N tracks: " << nTracks << std::endl;

    // Loop over tracker hits
    for (size_t itTrack = 0; itTrack < nTracks; itTrack++)
    {
        // Get the track
        EVENT::Track *track = static_cast<EVENT::Track *>(trackCollection->getElementAt(itTrack));

        // Loop over all hits in a track and mark them as used
        for (EVENT::TrackerHit *hit : track->getTrackerHits())
        {
            UsedHitsCollection->addElement(hit);
        }
    }

    // Add the hits to the output
    int nHits = trackerHitCollection->getNumberOfElements();
    // Set the map of responses 
    memset(&m_used, false, nHits);

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

    int nUsedHits = UsedHitsCollection->getNumberOfElements();

    streamlog_out(DEBUG4) << "  Total hits: " << nHits
                          << "  Used hits:  " << nUsedHits << std::endl;

    for (size_t itHit = 0; itHit < nUsedHits; itHit++)
    {
        // Get the hit
        TrackerHitPlane *hit = static_cast<TrackerHitPlane *>(UsedHitsCollection->getElementAt(itHit));
        unsigned int layer = myCellIDEncoding(hit)["layer"];
        unsigned int side = myCellIDEncoding(hit)["side"];
        unsigned int ladder = myCellIDEncoding(hit)["module"];
        unsigned int module = myCellIDEncoding(hit)["sensor"];

        streamlog_out(DEBUG0) << "Hit position " << layer << " " << ladder << " " << module << std::endl;
        const MySensorPos theOtherPos = {layer, side, ladder, module};

        // Checking if there are any hits in the other layer
        if (m_hitsMap.find(theOtherPos) == m_hitsMap.end())
        {
            streamlog_out(ERROR) << "No hit found, this should never happen!" << std::endl;
            continue;
        }

        for (size_t jitHit : m_hitsMap.at(theOtherPos))
        {
            TrackerHitPlane *hit2 = static_cast<TrackerHitPlane *>(trackerHitCollection->getElementAt(jitHit));
            TVector3 pos2(hit2->getPosition()[0], hit2->getPosition()[1], hit2->getPosition()[2]);

            if ((hit->getU() == hit2->getU()) && (hit->getV() == hit2->getV()))
            {
                streamlog_out(DEBUG0) << " --> found hit " << std::endl;
                m_used[itHit] = true;
            }
        }
    }

    // Once more to add the hits to the output (in slices)
    for (size_t itHit = 0; itHit < nHits; itHit++)
    {
        TrackerHitPlane *hit = static_cast<TrackerHitPlane *>(trackerHitCollection->getElementAt(itHit));

        if (!m_used[itHit])
        {
            SlimmedHitsCollection->addElement(hit);
        }
    }

    streamlog_out(DEBUG4) << "  Unused hits:  " << SlimmedHitsCollection->getNumberOfElements() << std::endl;

    // Store the filtered hit collections
    evt->addCollection(SlimmedHitsCollection, m_outputHitCollection);

    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    streamlog_out(DEBUG4) << "   done processing event: " << evt->getEventNumber()
                         << "   in run:  " << evt->getRunNumber() << std::endl;

    _nEvt++;
}

void HitSlimmer::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void HitSlimmer::end()
{

    //   std::cout << "HitSlimmer::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void HitSlimmer::getCollection(LCCollection *&collection, std::string collectionName, LCEvent *evt)
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
