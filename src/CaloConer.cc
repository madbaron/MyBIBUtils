#include "CaloConer.h"
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <filesystem>

#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <UTIL/LCRelationNavigator.h>
#include <IMPL/LCCollectionVec.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/CellIDEncoder.h>
#include <IMPL/LCRelationImpl.h>

#include "TLorentzVector.h"
#include "TVector3.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio;
using namespace marlin;

CaloConer aCaloConer;

CaloConer::CaloConer() : Processor("CaloConer")
{

    // Modify processor description
    _description = "CaloConer keeps only calo hits within fixed deltaR of MC truth";

    // Input collection
    registerProcessorParameter("MCParticleCollectionName",
                               "Name of the MCParticle input collection",
                               m_inputMCParticleCollection,
                               std::string("MCParticle"));

    // Input collection
    registerProcessorParameter("CaloHitCollectionName",
                               "Name of the CalorimeterHit input collection",
                               m_inputHitCollection,
                               std::string("EcalBarrelCollectionRec"));

    // Output collection
    registerProcessorParameter("GoodHitCollection",
                               "Good hits from calo",
                               m_outputHitCollection,
                               std::string("EcalBarrelCollectionConed"));

    // Input relation collection
    registerProcessorParameter("CaloRelationCollectionName",
                               "Name of the CalorimeterHit input relation collection",
                               m_inputRelationCollection,
                               std::string("EcalBarrelRelationsSimRec"));

    // Output relation collection
    registerProcessorParameter("GoodRelationCollection",
                               "Good hits SimRec relations",
                               m_outputRelationCollection,
                               std::string("EcalBarrelRelationsSimConed"));

    // Cone size
    registerProcessorParameter("ConeWidth",
                               "Cut in radians",
                               m_ConeSize,
                               0.2);
}

void CaloConer::init()
{

    streamlog_out(DEBUG) << "   init called  " << std::endl;

    // usually a good idea to
    printParameters();

    _nRun = 0;
    _nEvt = 0;
}

void CaloConer::processRunHeader(LCRunHeader *run)
{

    _nRun++;
}

void CaloConer::processEvent(LCEvent *evt)
{

    streamlog_out(DEBUG) << "Processing event " << _nEvt << std::endl;
    streamlog_out(DEBUG) << " in " << this->name() << std::endl;

    // Get the collection of MCParticles
    LCCollection *MCpartCollection = 0;
    getCollection(MCpartCollection, m_inputMCParticleCollection, evt);

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
            CalorimeterHit *hit = dynamic_cast<CalorimeterHit *>(caloHitCollection->getElementAt(itHit));

            // hit position
            TVector3 hitPos(hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);
            bool save = false;

            // Loop over MC particles and keep hit if within cone from particle 
            int nParts = MCpartCollection->getNumberOfElements();
            
            for (int itPart = 0; itPart < nParts; itPart++)
            {
                // Get MC part
                MCParticle *part = dynamic_cast<MCParticle *>(MCpartCollection->getElementAt(itPart));
                
                // --- Keep only the generator-level particles:
                if ( part->getGeneratorStatus() != 1 ) continue;
                
                TLorentzVector part_TLV(part->getMomentum()[0],part->getMomentum()[1],part->getMomentum()[2],part->getEnergy());

                double deltaR = fabs(part_TLV.Angle(hitPos));
                if(deltaR < m_ConeSize){
                    save =  true;
                    break;
                }
                
            }

            if (save)
            {
                streamlog_out(DEBUG0) << " accepted hit " << std::endl;

                outputHitCol->addElement(hit);

                LCRelation *rel = static_cast<LCRelation *>(inputHitRel->getElementAt(itHit));
                SimCalorimeterHit *simhit = static_cast<SimCalorimeterHit *>(rel->getTo());

                thitNav.addRelation(hit, simhit);
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

void CaloConer::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void CaloConer::end()
{
    //   std::cout << "CaloConer::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}

void CaloConer::getCollection(LCCollection *&collection, const std::string &collectionName, LCEvent *evt)
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
