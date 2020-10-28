import FWCore.ParameterSet.Config as cms

process = cms.Process("NtupleProduction")
isMC = True

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#from os import listdir
#DirPath = '/pnfs/knu.ac.kr/data/cms/store/user/sako/HAA4e/H800A3/HAA4e_H800A3_pythia8_AOD/HAA4e_H800A3_pythia8_AOD/180327_195007/0000/'
#InputFilesList = listdir(DirPath)
#InputFilesList = ['dcap://cluster142.knu.ac.kr/'+DirPath+x for x in InputFilesList if '.root' in x]

## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '' # InputFilesList
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string("80X_mcRun2_asymptotic_2016_TrancheIV_v6")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Ntuple.root')
)

process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) < 24 && position.Rho < 2"), # tracksSize() > 3 for the older cut
    filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)
# ElectronIDs = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']
# for idmod in ElectronIDs:
# setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# JEC
# process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
# process.ak4PFCorrectedJets = cms.EDProducer('CorrectedPFJetProducer',
#     src         = cms.InputTag('ak4PFJets'),
#     correctors  = cms.VInputTag('ak4PFL1L2L3Corrector')
# )

# newTask = cms.Task()
# process.egmGsfElectronIDSequence.associate(newTask)
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
process.load("Analysis.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi")
process.load("Analysis.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi")
# newTask.add(process.ModifiedHEEPIDVarValueMaps)
# newTask.add(process.ModifiedEcalRecHitIsolationScone)
# newTask.add(process.ModifiedHcalDepth1TowerIsolationScone)

process.tree = cms.EDAnalyzer("Ntuplizer",
    isMC = cms.untracked.bool(True),
    TriggerResults = cms.InputTag("TriggerResults","","HLT"),
    TriggerSummary = cms.InputTag("hltTriggerSummaryAOD"),
    Vertex = cms.InputTag("offlinePrimaryVertices"),
    PU = cms.InputTag("addPileupInfo"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    Muons = cms.InputTag("muons"),
    Electrons = cms.InputTag("gedGsfElectrons"),
    Conversions = cms.InputTag("conversions"),
    GenParticles = cms.InputTag("genParticles"),
    GenJets = cms.InputTag("ak4GenJets"),
    Generator = cms.InputTag("generator"),
    Photons = cms.InputTag("gedPhotons"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    miniRho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    # Jets = cms.InputTag("CorrectedPFJetProducer"),
    Jets = cms.InputTag("ak4PFJets"),
    MET = cms.InputTag("pfMet"),
    GsfTracks = cms.InputTag("electronGsfTracks"),
    GeneralTracks = cms.InputTag("generalTracks"),
    # HEEPVID = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
    nrSatCrysMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleNrSaturateIn5x5"),
    trkIsoMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleTrkPtIso"),
    nrMatchedTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleNrMatchedTrk"),
    rtMatchedTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleRtMatchedTrk"),
    addGsfTrkSelMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrkSel"),
    EcalRecHitIsoMap = cms.InputTag("ModifiedEcalRecHitIsolationScone","EcalRecHitIso"),
    noSelectedGsfTrk = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleNoSelectedGsfTrk"),
    addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"),
    KFParameters = cms.PSet(
        maxDistance = cms.double(0.01),
        maxNbrOfIterations = cms.int32(10)
    )
)

process.p = cms.Path(
    process.goodOfflinePrimaryVertices*
    # process.ak4PFL1L2L3CorrectorChain*
    # process.ak4PFCorrectedJets*
    # process.egmGsfElectronIDSequence*
    process.ModifiedHEEPIDVarValueMaps*
    process.ModifiedEcalRecHitIsolationScone*
    process.tree
)
