//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: IronFilterDetectorConstruction.hh $
//
/// \file IronFilterDetectorConstruction.hh
/// \brief Definition of the IronFilterDetectorConstruction class

#ifndef IronFilterDetectorConstruction_h
#define IronFilterDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class IronFilterDetectorMessenger;//pratyush

class IronFilterDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    IronFilterDetectorConstruction();
    virtual ~IronFilterDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();

    // get methods
    //
    ////const G4VPhysicalVolume* GetairlayersolidPV() const;
    ////const G4VPhysicalVolume* GetroomsolidPV() const;

    //lab
    const G4VPhysicalVolume* GetLabFloorWallsolidPV() const;
    const G4VPhysicalVolume* GetLabFloorExtendedsolidPV() const;
    const G4VPhysicalVolume* GetfrontglassdoorPV() const;
    const G4VPhysicalVolume* GetfrontdoorPV() const;
    const G4VPhysicalVolume* GetglasswindowPV() const;
    const G4VPhysicalVolume* GetreardoorPV() const;
    const G4VPhysicalVolume* GetInsulationPV() const;


    //fridge
    const G4VPhysicalVolume* GetDilutionUnitPV() const;
    const G4VPhysicalVolume* GetCPDUnitPV() const;
    const G4VPhysicalVolume* GetCopplerPillarCPDPV() const;
    const G4VPhysicalVolume* GetMixingPlatePV() const;
    const G4VPhysicalVolume* GetColdPlate2PV() const;
    const G4VPhysicalVolume* GetColdPlate1PV() const;
    const G4VPhysicalVolume* GetOneKPlatePV() const;
    const G4VPhysicalVolume* GetOneKShieldPV() const;
    const G4VPhysicalVolume* GetDilutionChamberPV() const;
    const G4VPhysicalVolume* GetFourKPlatePV() const;
    const G4VPhysicalVolume* GetFourKShieldPV() const;
    const G4VPhysicalVolume* GetSeventyKPlatePV() const;
    const G4VPhysicalVolume* GetSeventyKShieldPV() const;
    const G4VPhysicalVolume* GetOVCShieldPV() const;


  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    ////G4VPhysicalVolume* airlayer_solid_PV;
    ////G4VPhysicalVolume* room_solid_PV;

    //lab
    G4VPhysicalVolume* LabFloorWall_solid_PV;
    G4VPhysicalVolume* LabFloorExtended_solid_PV;
    G4VPhysicalVolume* frontglassdoor_PV;
    G4VPhysicalVolume* frontdoor_PV;
    G4VPhysicalVolume* glasswindow_PV;
    G4VPhysicalVolume* reardoor_PV;
    G4VPhysicalVolume* Insulation_PV;

    //fridge
    G4VPhysicalVolume* DilutionUnit_PV;
    G4VPhysicalVolume* CPD_Unit_PV;
    G4VPhysicalVolume* Coppler_Pillar_CPD_PV;
    G4VPhysicalVolume* MixingPlate_PV;
    G4VPhysicalVolume* ColdPlate2_PV;
    G4VPhysicalVolume* ColdPlate1_PV;
    G4VPhysicalVolume* OneKPlate_PV;
    G4VPhysicalVolume* OneKShield_PV;
    G4VPhysicalVolume* DilutionChamber_PV;
    G4VPhysicalVolume* FourKPlate_PV;
    G4VPhysicalVolume* FourKShield_PV;
    G4VPhysicalVolume* SeventyKPlate_PV;
    G4VPhysicalVolume* SeventyKShield_PV;
    G4VPhysicalVolume* OVCShield_PV;


    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
};



//////////********* LAB ***********//////////////
inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetLabFloorWallsolidPV() const {
  return LabFloorWall_solid_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetLabFloorExtendedsolidPV() const {
  return LabFloorExtended_solid_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetfrontglassdoorPV() const {
  return frontglassdoor_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetfrontdoorPV() const {
  return frontdoor_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetglasswindowPV() const {
  return glasswindow_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetreardoorPV() const {
  return reardoor_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetInsulationPV() const {
  return Insulation_PV;
}


////////******Fridge*********////////////////
inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetDilutionUnitPV() const {
  return DilutionUnit_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetCopplerPillarCPDPV() const {
  return Coppler_Pillar_CPD_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetCPDUnitPV() const {
  return CPD_Unit_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetMixingPlatePV() const {
  return MixingPlate_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetColdPlate2PV() const {
  return ColdPlate2_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetColdPlate1PV() const {
  return ColdPlate1_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetOneKPlatePV() const {
  return OneKPlate_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetOneKShieldPV() const {
  return OneKShield_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetDilutionChamberPV() const {
  return DilutionChamber_PV;
}


inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetFourKPlatePV() const {
  return FourKPlate_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetFourKShieldPV() const {
  return FourKShield_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetSeventyKPlatePV() const {
  return SeventyKPlate_PV;
}


inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetSeventyKShieldPV() const {
  return SeventyKShield_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetOVCShieldPV() const {
  return OVCShield_PV;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
