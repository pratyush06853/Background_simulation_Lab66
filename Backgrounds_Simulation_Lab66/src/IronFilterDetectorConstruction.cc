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
// $Id: IronFilterDetectorConstruction.cc $
//
/// \file IronFilterDetectorConstruction.cc
/// \brief Implementation of the IronFilterDetectorConstruction class

#include "IronFilterDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4NeutronHPThermalScatteringNames.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4AutoDelete.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Polycone.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDetectorConstruction::IronFilterDetectorConstruction()
 : G4VUserDetectorConstruction(),
   LabFloorWall_solid_PV(0),
   LabFloorExtended_solid_PV(0),
   frontglassdoor_PV(0),
   frontdoor_PV(0),
   glasswindow_PV(0),
   reardoor_PV(0),
   Insulation_PV(0),
   DilutionUnit_PV(0),
   CPD_Unit_PV(0),
   Coppler_Pillar_CPD_PV(0),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDetectorConstruction::~IronFilterDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* IronFilterDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IronFilterDetectorConstruction::DefineMaterials()
{
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density, fractionMass;
  G4String symbol, name;
  G4int nComponents, nAtoms;
  G4double temp;
  G4NistManager* NistMgr = G4NistManager::Instance();

  G4Element* elH  = new G4Element(name = "Hydrogen", symbol = "H", z = 1.0, a = 1.008*g/mole);
  G4Element* elLi6 = new G4Element(name = "Lithium6", symbol = "Li", z = 3.0, a = 6.015122*g/mole);
  G4Element* elLi7 = new G4Element(name = "Lithium6", symbol = "Li", z = 3.0, a = 7.016003*g/mole);
  G4Element* elB10 = new G4Element(name = "Boron10", symbol = "B", z = 5.0, a = 10.00*g/mole);
  G4Element* elB11 = new G4Element(name = "Boron11", symbol = "B", z = 5.0, a = 11.00*g/mole);
  G4Element* elB  = new G4Element(name = "Boron", symbol = "B", z = 5.0, a = 10.811*g/mole);
  G4Element* elC  = new G4Element(name = "Carbon", symbol = "C", z = 6.0, a = 12.011*g/mole);
  G4Element* elN = new G4Element("Nitrogen", symbol = "N",z = 7.,a = 14.01*g/mole);
  G4Element* elO  = new G4Element(name = "Oxygen", symbol = "O", z = 8.0, a = 15.999*g/mole);
  G4Element* elF = new G4Element(name= "Fluorine", symbol = "F", z = 9.0, a= 18.998403*g/mole); //pratyush
  G4Element* elNa = new G4Element("Sodium",symbol ="Na",z = 11.,a= 22.99*g/mole);
  G4Element* elMg = new G4Element("Magnesium",symbol ="Mg",z = 12.,a= 24.305*g/mole);
  G4Element* elAl  = new G4Element(name = "Aluminum", symbol = "Al", z = 13.0, a = 26.982*g/mole);
  G4Element* elSi  = new G4Element("Silicon", symbol = "Si", z = 14.0, a = 28.085*g/mole);
  G4Element* elP = new G4Element("Phosphorus",symbol ="P",z = 15.,a= 30.974*g/mole);
  G4Element* elS  = new G4Element(name = "sulphur", symbol = "S", z = 16.0, a = 32.065*g/mole);
  G4Element* elCl = new G4Element("Chlorine",symbol ="Cl",z = 17.,a= 35.453*g/mole);
  G4Element* elK = new G4Element("Potassium",symbol ="K",z = 19.,a= 39.098*g/mole);
  G4Element* elCa = new G4Element("Calcium",symbol ="Ca",z = 20.,a= 40.08*g/mole);
  G4Element* elFe  = new G4Element("Iron",symbol ="Fe",z = 26.,a= 55.85*g/mole);
  G4Element* elZn  = new G4Element(name = "zinc", symbol = "Ti", z = 30.0, a = 65.38*g/mole);
  G4Element* elRb = new G4Element("Rb",symbol ="Rb",z = 37.,a= 85.47 *g/mole);
  G4Element* elSr = new G4Element("Sr",symbol ="Sr",z = 38.,a= 87.62 *g/mole);
  G4Element* elZr = new G4Element("Zr",symbol ="Zr",z = 40.,a= 91.22 *g/mole);
  G4Element* elI  = new G4Element(name = "Iodine", symbol = "I", z = 53.0, a = 126.904473*g/mole);
  G4Element* elPb = new G4Element("Lead",symbol ="Pb", z = 82.,a= 207.19 *g/mole);

  G4Element *TS_H_P = new G4Element("TS_H_of_Polyethylene", "H", 1, 1.007*g/mole);
  G4Element *TS_H_W = new G4Element("TS_H_of_Water", "H", 1, 1.007*g/mole);


  //G4Element* elBe = new G4Element(name = "Beryllium", symbol = "Be", z = 4.0, a = 9.012*g/mole);
  //G4Element* elLi  = new G4Element(name = "Lithium", symbol = "Li", z = 3.0, a = 6.015*g/mole);
  //G4Element* elCr  = new G4Element(name = "Chromium", symbol = "Cr", z = 24.0, a = 51.996*g/mole);
  //G4Element* elFe  = new G4Element(name = "Iron", symbol = "Fe", z = 26.0, a = 55.845*g/mole);
  //G4Element* elNi  = new G4Element(name = "Nickel", symbol = "Ni", z = 28.0, a = 58.693*g/mole);
  //G4Element* elMo  = new G4Element(name = "Molybdenum", symbol = "Mo", z = 42.0, a = 95.94*g/mole);
  //G4Element* elI  = new G4Element(name = "Iodine", symbol = "I", z = 53.0, a = 127*g/mole);


  //Vacuum
  new G4Material("galactic", z = 1.0, a = 1.01*g/mole, density = universe_mean_density, kStateGas, 2.73*kelvin, 3.e-18*pascal);


  //helium
  new G4Material("liquid_helium", z = 2., a = 4.00*g/mole, density = 0.145*g/cm3, kStateLiquid, temp = 3.*kelvin);

  // Lead
  new G4Material(name = "Pb", z = 82.0, a = 207.2*g/mole, density = 11.34*g/cm3);
  new G4Material("NatLead", 82.0, 207.19*g/mole,  11.36*g/cm3, kStateSolid, 296*kelvin);

  //Copper
  new G4Material("NatCu", z =29.,  a =63.55*g/mole,density =   8.96 *g/cm3);

  // Iron
  new G4Material("NatIron", z = 26.0, a = 55.845*g/mole, density = 7.874*g/cm3);

  // scandium
  new G4Material("NatScandium", z = 21.0, a = 44.95*g/mole, density = 2.985*g/cm3);

  //Aluminum
  new G4Material("NatAluminum", z = 13.0, a = 26.9815384*g/mole, density = 2.70*g/cm3, kStateSolid, 296*kelvin);


  //Natural B
  G4Element* NatB=NistMgr->FindOrBuildElement("B");
  //Natural C
  G4Element* NatC=NistMgr->FindOrBuildElement("C");
  //Natural H
  G4Element* NatH=NistMgr->FindOrBuildElement("H");
  //Natural Oxygen
  G4Element* NatO=NistMgr->FindOrBuildElement("O");
  //Natural Na
  G4Element* NatNa=NistMgr->FindOrBuildElement("Na");


  // Silicon
  G4Isotope* Si28 = new G4Isotope("Si28", 14, 28, 27.9769265350*g/mole);
  G4Isotope* Si29 = new G4Isotope("Si29", 14, 29, 28.9764946653*g/mole);
  G4Isotope* Si30 = new G4Isotope("Si30", 14, 30, 29.973770137*g/mole);

  G4Element* NatSi = new G4Element("silicon_natural", "Si",3);
  NatSi->AddIsotope(Si28,92.2*perCent);
  NatSi->AddIsotope(Si29,4.7*perCent);
  NatSi->AddIsotope(Si30,3.1*perCent);

  G4Material* silicon = new G4Material("silicon", 2.3290*g/cm3,1, kStateSolid, 296*kelvin);
  silicon->AddElement(NatSi, 1);


  // Titanium
  G4Isotope* Ti46 = new G4Isotope("Ti46", 22, 46, 45.9526316*g/mole);
  G4Isotope* Ti47 = new G4Isotope("Ti47", 22, 47, 46.9517631*g/mole);
  G4Isotope* Ti48 = new G4Isotope("Ti48", 22, 48, 47.9479463*g/mole);
  G4Isotope* Ti49 = new G4Isotope("Ti49", 22, 49, 48.9478700*g/mole);
  G4Isotope* Ti50 = new G4Isotope("Ti50", 22, 50, 49.9447912*g/mole);

  G4Element* NatTi = new G4Element("titanium_natural", "Ti",5);
  NatTi->AddIsotope(Ti46,8.25*perCent);
  NatTi->AddIsotope(Ti47,7.44*perCent);
  NatTi->AddIsotope(Ti48,73.72*perCent);
  NatTi->AddIsotope(Ti49,5.41*perCent);
  NatTi->AddIsotope(Ti50,5.18*perCent);

  G4Material* titanium = new G4Material("titanium", 4.507*g/cm3,1, kStateSolid, 296*kelvin);
  titanium->AddElement(NatTi, 1);



  //Ni60
  G4Material* Ni60 = new G4Material("Ni60", z = 28.0, a = 59.9307864*g/mole, density = 8.9*g/cm3);
  //Fe54
  G4Material* Fe54 = new G4Material("Fe54", z = 26.0, a = 	53.9396090*g/mole, density = 7.874*g/cm3);
  //Co
  G4Material* NatCo = new G4Material("NatCo", z = 27.0, a = 	58.933*g/mole, density = 8.86*g/cm3);





  //LiI
  G4Isotope* Li6 = new G4Isotope("Li6", 3, 6, 6.01*g/mole);
  G4Isotope* Li7 = new G4Isotope("Li7", 3, 7, 7.01*g/mole);

  G4Element* Li6enriched = new G4Element("Lithium_enchried", "Li",2);
  Li6enriched->AddIsotope(Li6,97.0*perCent);
  Li6enriched->AddIsotope(Li7,3.0*perCent);

  G4Material* Li6F = new G4Material("Li6F", density= 2.64 * g / cm3,nComponents= 2, kStateSolid, 296*kelvin);
  Li6F->AddElement(Li6enriched, 1);
  Li6F->AddElement(elF,1);

  G4Material* ZnS = new G4Material("ZnS", density= 4.09 * g / cm3,nComponents= 2, kStateSolid, 296*kelvin);
  ZnS->AddElement(elZn, 1);
  ZnS->AddElement(elS,1);

  //Remeber HD has 1:2 LiF to ZnS concentration
  G4Material* ej_426_HD= new G4Material( "ej_426_HD", density=3.60*g/cm3, nComponents=2, kStateSolid, 296*kelvin);
  ej_426_HD->AddMaterial( Li6F, 33.3*perCent );   //1.026*g/cm3
  ej_426_HD->AddMaterial( ZnS, 66.7*perCent );


  //Polyethylene moderator
  G4Material*  polyethylene = new G4Material("polyethylene", density=0.94*g/cm3, nComponents=2,kStateSolid, 296*kelvin);
  polyethylene->AddElement(NatC, 1);
  polyethylene->AddElement(TS_H_P, 2);

  // Assuming PMMA -- see
  //	http://en.wikipedia.org/wiki/Poly(methyl_methacrylate)
  G4Material*  acrylic = new G4Material("acrylic", density= 1.17 * g/cm3, nComponents=3 ,kStateSolid, 296*kelvin);
  acrylic->AddElement(NatC, 5);
  acrylic->AddElement(elO, 2);
  acrylic->AddElement(TS_H_P, 8);

  //Water
  G4Material* water = new G4Material("water", density= 1.00 * g / cm3,nComponents= 2, kStateLiquid, 296*kelvin);
  water->AddElement(TS_H_W,2);
  water->AddElement(elO,1);

  //Aluminum Fluoride
  G4Material* AlF3 = new G4Material("AlF3", density= 3.10 * g / cm3,nComponents= 2); //pratyush
  AlF3->AddElement(elAl, 1);  //pratyush
  AlF3->AddElement(elF, 3);   //pratyush

  //Borax
  G4Material* Borax = new G4Material("Borax", density= 0.76* g / cm3,nComponents= 4,kStateSolid, 296*kelvin); //pratyush
  Borax->AddElement(NatNa,12.06*perCent);//2
  Borax->AddElement(NatB,11.34*perCent);//4
  Borax->AddElement(NatH,5.29*perCent);//20
  Borax->AddElement(NatO,71.32*perCent);//17

  //Boric Acid https://www.convertunits.com/molarmass/Boric+Acid
  G4Material* boric_acid = new G4Material("boric_acid", density= 1.44* g / cm3,nComponents= 3,kStateSolid, 296*kelvin); //pratyush
  boric_acid->AddElement(NatH,4.890*perCent);//3
  boric_acid->AddElement(NatB,17.484*perCent);//1
  boric_acid->AddElement(NatO,77.626*perCent);//3

  //Borax_water_Mixture(5.8% solubity of borax https://omsi.edu/sites/all/FTP/files/kids/Borax-msds.pdf)
  //mixture of 5.5% Borax and 94.5% of Water
  G4Material* borax_water = new G4Material( "borax_water",density= 0.9868*g/cm3, nComponents= 2,kStateLiquid, 296*kelvin); //pratyush
  borax_water->AddMaterial( Borax, 5.5*perCent );
  borax_water->AddMaterial( water, 94.5*perCent );


  //Borax_BoricAcid_buffer(https://www.researchgate.net/publication/244069630_Preparation_of_highly_concentrated_aqueous_solution_of_sodium_borate)
  //mixture of 20g BoricAcid, 25g of Borax and 100g of water
  G4Material* borax_boricacid_buffer = new G4Material( "borax_boricacid_buffer",density= 1.019*g/cm3, nComponents= 3,kStateLiquid, 296*kelvin);
  borax_boricacid_buffer->AddMaterial( boric_acid, 13.7*perCent );//pratyush
  borax_boricacid_buffer->AddMaterial( Borax, 17.2*perCent );//pratyush
  borax_boricacid_buffer->AddMaterial( water, 69.1*perCent );//pratyush

  //Fluental
  //mixture of 40% Al and 60% of AlF_3
  G4Material* fluental = new G4Material( "fluental",density= 2.94*g/cm3, nComponents= 2); //pratyush
  fluental->AddMaterial( AlF3, 60.*perCent );
  fluental->AddElement( elAl, fractionMass = 40.*perCent );

  //polyethyleneBorated 3 % borated poly
  G4Material* boratedPoly = new G4Material( "boratedPoly", density=1.19*g/cm3, nComponents=3,kStateSolid, 296*kelvin);
  boratedPoly->AddElement( NatB, 3.*perCent );
  boratedPoly->AddElement( NatC, 82.576*perCent );
  boratedPoly->AddElement(TS_H_P, 14.424*perCent );

  //wood
  G4Material* wood = new G4Material("wood", density=0.9*g/cm3, nComponents=3);
  wood->AddElement(TS_H_P , 4);
  wood->AddElement(elO , 1);
  wood->AddElement(elC, 2);

  G4Material* quartz = new G4Material("quartz", density=2.200*g/cm3, nComponents=2);
  quartz->AddElement(elSi, 1);
  quartz->AddElement(elO , 2);

  G4Material*soft_tissue = new G4Material("soft_tissue",density= 0.9869*g/cm3,nComponents=9);
  soft_tissue->AddElement(elH,0.105);
  soft_tissue->AddElement(elC,0.256);
  soft_tissue->AddElement(elN,0.027);
  soft_tissue->AddElement(elO,0.602);
  soft_tissue->AddElement(elNa,0.001);
  soft_tissue->AddElement(elP,0.002);
  soft_tissue->AddElement(elS,0.003);
  soft_tissue->AddElement(elCl,0.002);
  soft_tissue->AddElement(elK,0.002);


  //soil
  G4Material*soil = new G4Material("soil",density= 1.50*g/cm3,nComponents=8);
  soil->AddElement(elH,0.021);
  soil->AddElement(elC,0.016);
  soil->AddElement(elO,0.577);
  soil->AddElement(elAl,0.050);
  soil->AddElement(elSi,0.271);
  soil->AddElement(elK,0.013);
  soil->AddElement(elCa,0.041);
  soil->AddElement(elFe,0.011);

  //concrete
  G4Material*concrete = new G4Material("concrete",density= 2.3*g/cm3,nComponents=10);
  concrete->AddElement(elH,0.01);
  concrete->AddElement(elC,0.001);
  concrete->AddElement(elO,0.529107);
  concrete->AddElement(elNa,0.016);
  concrete->AddElement(elMg,0.002);
  concrete->AddElement(elAl,0.033872);
  concrete->AddElement(elSi,0.337021);
  concrete->AddElement(elK,0.013);
  concrete->AddElement(elCa,0.044);
  concrete->AddElement(elFe,0.014);


  //borated concrete
  G4Material* borated_concrete = new G4Material("borated_concrete",density= 2.32*g/cm3,nComponents=10);
  borated_concrete->AddElement(elH,0.96*perCent);
  borated_concrete->AddElement(elC,5.36*perCent);
  borated_concrete->AddElement(elO,51.3*perCent);
  borated_concrete->AddElement(NatB,2.9*perCent);
  borated_concrete->AddElement(elMg,0.42*perCent);
  borated_concrete->AddElement(elAl,0.79*perCent);
  borated_concrete->AddElement(elSi,15.7*perCent);
  borated_concrete->AddElement(elS,0.42*perCent);
  borated_concrete->AddElement(elCa,23*perCent);
  borated_concrete->AddElement(elFe,0.50*perCent);



  // Print materials
  //G4cout <<"Is hydrogen declared in the scintillator a Thermal Hydrogen"  <<G4NeutronHPThermalScatteringNames::IsThisThermalElement(TS_H_of_Polyethylene)<< G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* IronFilterDetectorConstruction::DefineVolumes()
{

  // Get materials
  G4Material* Vacuum = G4Material::GetMaterial("galactic");
  G4Material* Helium = G4Material::GetMaterial("liquid_helium");
  G4Material* Titanium = G4Material::GetMaterial("titanium");
  G4Material* Aluminum = G4Material::GetMaterial("NatAluminum");
  G4Material* Lead = G4Material::GetMaterial("NatLead");
  G4Material* Silicon = G4Material::GetMaterial("silicon");
  G4Material* Iron = G4Material::GetMaterial("NatIron");
  G4Material* Scandium = G4Material::GetMaterial("NatScandium");
  G4Material* Fluental = G4Material::GetMaterial("fluental");
  G4Material* BoraxWater = G4Material::GetMaterial("borax_water");
  G4Material* BoraxBoricAcidBuffer = G4Material::GetMaterial("borax_boricacid_buffer");
  G4Material* Copper= G4Material::GetMaterial("NatCu");

  G4Material* EJ4265HD = G4Material::GetMaterial("ej_426_HD");
  G4Material* Polyethylene = G4Material::GetMaterial("polyethylene");
  G4Material* Acrylic = G4Material::GetMaterial("acrylic");

  G4Material*  Soft_Tissue=G4Material::GetMaterial("soft_tissue");
  G4Material*  Concrete = G4Material::GetMaterial("concrete");
  G4Material*  Borated_Concrete = G4Material::GetMaterial("borated_concrete");
  G4Material*  Wood = G4Material::GetMaterial("wood");
  G4Material*  Quartz = G4Material::GetMaterial("quartz");
  G4Material*  Soil = G4Material::GetMaterial("soil");


  if ( ! Vacuum ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("IronFilterDetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }



//
// Sizes
//




  G4double DD_Height = 20.0*cm;
  G4double DD_Extra_Height= 60.0*cm;
  G4double shieldthickness = 40.0*cm; //20.0*cm; //20
  G4double leadshieldthickness = 20.0*cm;
  G4double Polyshieldthickness = 40.0*cm;
  G4double SourceRadius = 1.7*cm;
  G4double delta= 1.0*cm;// 1.0cm parameter of catchment area



  ///////////***********lab room************//////////////////
  G4double colimator_length=26.0*cm;
  G4double fFilterCellSpacing= 50.0*cm+26.0*cm;

  G4double Water_cylindercal_can_radius = 152.7175*cm;
  G4double Water_cylindercal_can_height = 115.8875*cm;
  G4double ConcreteSupport_height = 80.0*cm;
  G4double DT_Ti_T_location = 207.5*mm;
  G4double Insulation_Thickness = 5*mm;


  G4double lab68_wall_thickness = 25.0*cm ;
  //distance from the outer egdes
  G4double lab68_wall_x = 6.654*m ;
  G4double lab68_wall_y = 9.35063*m ;
  G4double lab68_wall_z = 5.636*m ;

  G4double Pump_chase_y = 2.9*m ;

  G4double lab68_frontdoor_glass_height = 2.9*m;
  G4double lab68_frontdoor_glass_width = 2.2*m;

  G4double lab68_frontdoor_wood_height = 2.3*m;
  G4double lab68_frontdoor_wood_width = 1.5*m;

  G4double lab68_glasswindow_height = 1.08*m;
  G4double lab68_glasswindow_width = 0.57*m;

  G4double lab68_reardoor_height = 2.3*m;
  G4double lab68_reardoor_width = 0.91*m;

  G4double lab68_frontdoor_x_coordinate = lab68_wall_x/2 -lab68_frontdoor_glass_width/2.0-1.2*m;
  G4double lab68_reardoor_x_coordinate = -lab68_wall_x/2.0 +lab68_reardoor_width/2.0+2.613*m;

  G4ThreeVector position_of_origin = {2.721*m, -2.703*m, 1.4547*m}; //with repect to the outer upper left corner of the room(Doug's corner)

  G4ThreeVector xyposition_of_origin = {2.721*m, -2.703*m, 0};

  //test box that surrounds the fridge
  //G4double test_x= 50*cm;
  //G4double test_y= 50*cm;
  //G4double test_z= 2.3*m;


  //test box that surrounds the mixing plate
  G4double test_x= 40*cm;
  G4double test_y= 40*cm;
  G4double test_z= 10*cm;
  //Fridge
  G4double SeventyKShield_Width = 0.2*cm;
  G4double SeventyKShield_Radius = 18.95*cm;
  G4double SeventyKShield_Height = 1.28*m;

  G4double OVCShield_Width = 0.4*cm;
  G4double OVCShield_Radius = 21.1*cm;
  G4double OVCShield_Height = 1.5*m;

  G4double FourKShield_Width = 0.2*cm;
  G4double FourKShield_Radius = 17.2*cm;
  G4double FourKShield_Height = 1.048*m;

  G4double OneKShield_Width = 0.2*cm;
  G4double OneKShield_Radius = 15.65*cm;
  G4double OneKShield_Height = 0.881*m;

  //G4double DilutionUnit_Radius = 10.0*cm;
  G4double DilutionUnit_Radius = 3.0*cm;
  G4double DilutionUnit_Height = 2.75*cm;

  G4double DilutionChamber_Radius = DilutionUnit_Radius;//10.0*cm
  G4double DilutionChamber_Height = DilutionUnit_Height + 5.3*cm;//15.3*cm
  G4double DilutionChamber_Width = 2.0*mm;
  G4double DilutionChamber_bottomplate_thick = 5*mm;
  G4double DilutionChamber_upperplate_thick = 3*mm;

  G4double SeventyKPlate_Radius = 20.0*cm;
  G4double FourKPlate_Radius = 18.3*cm;
  G4double OneKPlate_Radius = 16.7*cm;
  G4double ColdPlate1_Radius = 15.0*cm;
  G4double ColdPlate2_Radius = 15.0*cm;
  G4double MixingPlate_Radius = 15.0*cm;

  G4double SeventyKPlate_Thickness = 0.4*cm;
  G4double FourKPlate_Thickness  = 0.4*cm;
  G4double OneKPlate_Thickness  = 0.4*cm;
  G4double ColdPlate1_Thickness  = 0.4*cm;
  G4double ColdPlate2_Thickness  = 0.4*cm;
  G4double MixingPlate_Thickness  = 0.5*cm;

  G4double OVC_z = 10.3*cm;
  G4double SeventyKPlate_z = 20.3*cm+ FourKPlate_Thickness/2.0 + SeventyKPlate_Thickness/2.0;
  G4double FourKPlate_z  = 14.1*cm+ OneKPlate_Thickness/2.0 + FourKPlate_Thickness/2.0;
  G4double OneKPlate_z  = 9.4*cm + ColdPlate1_Thickness/2.0 + OneKPlate_Thickness/2.0;
  G4double ColdPlate1_z  = 9.6*cm + ColdPlate1_Thickness/2.0+ ColdPlate2_Thickness/2.0;
  G4double ColdPlate2_z  = 16*cm + MixingPlate_Thickness/2.0+ ColdPlate2_Thickness/2.0;
  //G4double MixingPlate_z  = 10.4*cm + MixingPlate_Thickness/2;
  G4double MixingPlate_z  = 20.7*cm + MixingPlate_Thickness/2;


  ///backing Detector
  G4double zeroRadius = 0.*cm;
  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;
  //G4double shieldHeight = 400.*cm;



//
// Rotations
//

  G4RotationMatrix* NO_ROT = new G4RotationMatrix;
  G4RotationMatrix* xRot7 = new G4RotationMatrix;
  G4RotationMatrix* xRot8 = new G4RotationMatrix;
  xRot7 -> rotateX(-pi/2.*rad);
  xRot8 -> rotateX(pi/2.*rad);

  G4RotationMatrix* turnAlongX = new G4RotationMatrix;
  turnAlongX->rotateX(90*deg);

  G4RotationMatrix* turnAlongXY = new G4RotationMatrix;
  turnAlongXY->rotateZ(90*deg);
  turnAlongXY->rotateZ(-100*deg);

  G4RotationMatrix* turnAlong = new G4RotationMatrix;
  turnAlong->rotateZ(10*deg);
  //turnAlong->rotateY(10*deg);

  G4RotationMatrix* turnAlong190 = new G4RotationMatrix;
  turnAlong190->rotateZ(190*deg);

  G4RotationMatrix* turnAlongZ = new G4RotationMatrix;
  turnAlongZ->rotateZ(90*deg);
  turnAlongZ->rotateZ(-90*deg);

//
// GEOMETRY
//




  //G4VSolid* vacuum_solid = new G4Tubs("vacuum_solid", zeroRadius,Inner_Radius+ Radial_thickness+20.0*cm, 100.0*cm, startAngle, spanningAngle);
  G4VSolid* vacuum_solid =new G4Box("vacuum_solid", 150.0*m, 150.0*m, 150.0*m);
  //G4VSolid* vacuum_solid = new G4Box("vacuum_solid", CapturerLength/2.0, CapturerLength/2.0, (shieldHeight)/2.0);
  G4LogicalVolume* vacuum_solid_LV = new G4LogicalVolume(vacuum_solid, Vacuum, "vacuum_solid");
  G4VPhysicalVolume* vacuum_solid_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(), vacuum_solid_LV, "Vacuum_solid", 0, false, 0, fCheckOverlaps);


  //Lab66
  G4VSolid* Main_2_S = new G4Box("Main_2_solid", lab68_wall_x/2.0, lab68_wall_y/2.0 , lab68_wall_z/2.0);
  G4VSolid* hole_2_S = new G4Box("hole_2_solid", (lab68_wall_x-2*lab68_wall_thickness)/2.0, (lab68_wall_y-2*lab68_wall_thickness)/2.0, (lab68_wall_z-2*lab68_wall_thickness)/2.0);
  G4SubtractionSolid* LabFloorWall_solid_S= new G4SubtractionSolid("LabFloorWall_solid", Main_2_S, hole_2_S, NO_ROT, G4ThreeVector(0.,0., 0.));
  G4LogicalVolume* LabFloorWall_solid_LV = new G4LogicalVolume(LabFloorWall_solid_S, Concrete, "LabFloorWall_solid");
  LabFloorWall_solid_PV = new G4PVPlacement(turnAlong, G4ThreeVector{lab68_wall_x/2.0,-lab68_wall_y/2.0,lab68_wall_z/2.0}-position_of_origin, LabFloorWall_solid_LV, "LabFloorWall", vacuum_solid_LV, false, 0, fCheckOverlaps);
  LabFloorWall_solid_LV->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));

  G4VSolid* frontglassdoor_S = new G4Box("frontglassdoor_solid", lab68_frontdoor_glass_width/2.0, lab68_wall_thickness/2.0, lab68_frontdoor_glass_height/2.0);
  G4LogicalVolume* frontglassdoor_LV = new G4LogicalVolume(frontglassdoor_S, Quartz, "frontglassdoor");
  frontglassdoor_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(lab68_frontdoor_x_coordinate ,-(lab68_wall_y-lab68_wall_thickness)/2.0, -(lab68_wall_z-2*lab68_wall_thickness-lab68_frontdoor_glass_height)/2.0), frontglassdoor_LV, "Front Glass Door", LabFloorWall_solid_LV, false, 0, fCheckOverlaps);
  frontglassdoor_LV->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));


  G4VSolid* frontdoor_S = new G4Box("frontdoor_solid", lab68_frontdoor_wood_width/2.0, lab68_wall_thickness/2.0, lab68_frontdoor_wood_height/2.0);
  //G4SubtractionSolid* Main_2b_S= new G4SubtractionSolid("Main_2b_solid", Main_2a_S, frontdoor_S, NO_ROT, G4ThreeVector(lab68_frontdoor_x_coordinate ,-lab68_wall_y/2.0, -(lab68_wall_z-2*lab68_wall_thickness-lab68_frontdoor_glass_height)/2.0));
  G4LogicalVolume* frontdoor_LV = new G4LogicalVolume(frontdoor_S, Wood, "frontdoor");
  frontdoor_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(lab68_frontdoor_glass_width/2.0-lab68_frontdoor_wood_width/2.0 ,0, -lab68_frontdoor_glass_height/2.0+lab68_frontdoor_wood_height/2.0), frontdoor_LV, "Front Wood Door", frontglassdoor_LV, false, 0, fCheckOverlaps);
  frontdoor_LV->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));


  G4VSolid* glasswindow_S = new G4Box("glasswindow_solid", lab68_glasswindow_width/2.0, lab68_wall_thickness/2.0, lab68_glasswindow_height/2.0);
  G4LogicalVolume* glasswindow_LV = new G4LogicalVolume(glasswindow_S, Quartz, "glasswindow");
  glasswindow_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(lab68_frontdoor_wood_width/2.0-0.16*m-lab68_glasswindow_width/2.0, 0, -lab68_frontdoor_wood_height/2.0+1.11*m+lab68_glasswindow_height/2.0), glasswindow_LV, "Front Glass Window", frontdoor_LV, false, 0, fCheckOverlaps);
  glasswindow_LV->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));


  G4VSolid* reardoor_S = new G4Box("reardoor_solid", lab68_reardoor_width/2.0, lab68_wall_thickness/2.0, lab68_reardoor_height/2.0);
  G4LogicalVolume* reardoor_LV = new G4LogicalVolume(reardoor_S, Wood, "reardoor");
  reardoor_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(lab68_reardoor_x_coordinate ,(lab68_wall_y-lab68_wall_thickness)/2.0, -(lab68_wall_z-2*lab68_wall_thickness-lab68_reardoor_height)/2.0), reardoor_LV, "Rear Wooden Door", LabFloorWall_solid_LV, false, 0, fCheckOverlaps);
  reardoor_LV->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));

  G4double soil_width=2*m;
  //G4VSolid* LabFloorExtended_solid_S=  new G4Box("LabFloorExtended_solid", 20.0*m, 20.0*m , soil_width/2.0);
  G4VSolid* LabFloorExtended_solid_S=  new G4Box("LabFloorExtended_solid", lab68_wall_x/2.0, lab68_wall_y/2.0  , soil_width/2.0);
  //G4LogicalVolume* LabFloorExtended_solid_LV = new G4LogicalVolume(LabFloorExtended_solid_S, Soil, "LabFloorExtended_solid");
  G4LogicalVolume* LabFloorExtended_solid_LV = new G4LogicalVolume(LabFloorExtended_solid_S, Concrete, "LabFloorExtended_solid");
  LabFloorExtended_solid_PV = new G4PVPlacement(turnAlong, G4ThreeVector{lab68_wall_x/2.0,-lab68_wall_y/2.0,lab68_wall_z/2.0}-position_of_origin-G4ThreeVector(0., 0., lab68_wall_z/2.0+soil_width/2.0), LabFloorExtended_solid_LV, "LabFloor_extended", vacuum_solid_LV, false, 0, fCheckOverlaps);
  LabFloorExtended_solid_LV->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));
  //LabFloorExtended_solid_LV->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));




///Lab Ends///////////
////////////////////////////******************///////////////////////////////////////////


////Fridge_Begins///////
///////////////////////////*******************///////////////////////////////////////////
//DilutionUnit contains superfluid Helium
G4VSolid* DilutionUnit_S = new G4Tubs( "DilutionUnit", zeroRadius, DilutionUnit_Radius, (DilutionUnit_Height /2.0), startAngle, spanningAngle);
G4LogicalVolume *DilutionUnit_LV = new G4LogicalVolume( DilutionUnit_S, Helium, "DilutionUnit" );
//G4LogicalVolume *DilutionUnit_LV = new G4LogicalVolume( DilutionUnit_S, Germanium, "DilutionUnit" );
DilutionUnit_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0,0,0),DilutionUnit_LV, "DilutionUnit", vacuum_solid_LV, false, 0, fCheckOverlaps);
DilutionUnit_LV->SetVisAttributes(G4VisAttributes(G4Colour::Green()));


//Dilution Chamber that holds the superfluid Helium
//G4double nedges[7]=  {-DilutionUnit_Height /2.0-DilutionChamber_bottomplate_thick,
//                         -(DilutionUnit_Height /2.0),-(DilutionUnit_Height /2.0),
//                         0.0,
//                         DilutionChamber_Height-DilutionUnit_Height/2.0, DilutionChamber_Height-DilutionUnit_Height/2.0,
//                         DilutionChamber_Height-DilutionUnit_Height/2.0+DilutionChamber_upperplate_thick};
//G4double innerradius[7]= {0.0, 0.0, DilutionChamber_Radius ,DilutionChamber_Radius ,  DilutionChamber_Radius, 0.0, 0.0 };
//G4double outerradius[7]= {DilutionChamber_Radius+DilutionChamber_Width, DilutionChamber_Radius+DilutionChamber_Width, DilutionChamber_Radius+DilutionChamber_Width,
//  DilutionChamber_Radius+DilutionChamber_Width, DilutionChamber_Radius+DilutionChamber_Width, DilutionChamber_Radius+DilutionChamber_Width
//, DilutionChamber_Radius+DilutionChamber_Width};


//Dilution Chamber
G4double nedges[15]= { -48.75*mm,
                       -13.75*mm,
                       -13.75*mm,
                       0*mm,
                       13.75*mm,
                       13.75*mm,
                       38.75*mm,
                       38.75*mm,
                       52.25*mm,
                       52.25*mm,
                       64.75*mm,
                       64.75*mm,
                       104.25*mm,
                       104.25*mm,
                       107.25*mm };



  G4double innerradius[15]= { 0*mm,
                      0*mm,
                      30*mm,
                      30*mm,
                      30*mm,
                      52.5*mm,
                      52.5*mm,
                      18*mm,
                      18*mm,
                      11*mm,
                      11*mm,
                      98*mm,
                      98*mm,
                      0*mm,
                      0*mm};


      G4double outerradius[15]= { 100*mm,
              100*mm,
              100*mm,
              100*mm,
              100*mm,
              100*mm,
              100*mm,
              100*mm,
              100*mm,
              100*mm,
              100*mm,
             100*mm,
            100*mm,
            100*mm,
            100*mm};

G4ThreeVector position_DilutionChamber = G4ThreeVector(0, 0, 0);

G4Polycone* DilutionChamber_S = new G4Polycone("DilutionChamber", startAngle, spanningAngle, 15, nedges, innerradius, outerradius);
G4LogicalVolume*  DilutionChamber_LV= new G4LogicalVolume(DilutionChamber_S, Copper, "DilutionChamber");
DilutionChamber_PV = new G4PVPlacement(NO_ROT, position_DilutionChamber , DilutionChamber_LV, "DilutionChamber", vacuum_solid_LV, false, 0, fCheckOverlaps);
DilutionChamber_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));


G4VSolid* CPD_Unit_S = new G4Tubs( "CPD_Unit", zeroRadius, (76*mm/2.0), (1*mm /2.0), startAngle, spanningAngle);
G4LogicalVolume *CPD_Unit_LV = new G4LogicalVolume( CPD_Unit_S, Silicon, "CPD_Unit" );
//new G4PVPlacement( NO_ROT, G4ThreeVector(0,0,20.25*mm),CPD_Unit_LV, "CPD", vacuum_solid_LV, false, 0, fCheckOverlaps);
CPD_Unit_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0,0,19.75*mm),CPD_Unit_LV, "CPD", vacuum_solid_LV, false, 0, fCheckOverlaps);
CPD_Unit_LV->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));


G4VSolid* Coppler_Pillar_CPD_S = new G4Tubs( "Coppler_Pillar_CPD", zeroRadius, (17.8*mm/2.0), (83.5*mm /2.0), startAngle, spanningAngle);
G4LogicalVolume *Coppler_Pillar_CPD_LV = new G4LogicalVolume( Coppler_Pillar_CPD_S, Copper, "Coppler_Pillar_CPD" );
//new G4PVPlacement( NO_ROT, G4ThreeVector(0,0,20.75*mm+ 83.5*mm /2.0),Coppler_Pillar_CPD_LV, "CPD", vacuum_solid_LV, false, 0, fCheckOverlaps);
Coppler_Pillar_CPD_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0,0,20.75*mm+ 83.5*mm /2.0),Coppler_Pillar_CPD_LV, "Copper_Pillar", vacuum_solid_LV, false, 0, fCheckOverlaps);
Coppler_Pillar_CPD_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));



//MixingPlate
G4VSolid* MixingPlate_S = new G4Tubs( "MixingPlate", zeroRadius, MixingPlate_Radius, (MixingPlate_Thickness /2.0), startAngle, spanningAngle);
G4LogicalVolume *MixingPlate_LV = new G4LogicalVolume( MixingPlate_S, Copper, "MixingPlate" );
MixingPlate_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0, 0, MixingPlate_z),MixingPlate_LV, "MixingPlate", vacuum_solid_LV, false, 0, fCheckOverlaps);
MixingPlate_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));

//ColdPlate2
G4VSolid* ColdPlate2_S = new G4Tubs( "ColdPlate2", zeroRadius, ColdPlate2_Radius, (ColdPlate2_Thickness /2.0), startAngle, spanningAngle);
G4LogicalVolume *ColdPlate2_LV = new G4LogicalVolume( ColdPlate2_S, Copper, "ColdPlate2" );
ColdPlate2_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0, 0, MixingPlate_z+ColdPlate2_z),ColdPlate2_LV, "ColdPlate2", vacuum_solid_LV, false, 0, fCheckOverlaps);
ColdPlate2_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));

//ColdPlate1
G4VSolid* ColdPlate1_S = new G4Tubs( "ColdPlate1", zeroRadius, ColdPlate1_Radius, (ColdPlate1_Thickness /2.0), startAngle, spanningAngle);
G4LogicalVolume *ColdPlate1_LV = new G4LogicalVolume( ColdPlate1_S, Copper, "ColdPlate1" );
ColdPlate1_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0, 0, MixingPlate_z+ColdPlate2_z+ColdPlate1_z),ColdPlate1_LV, "ColdPlate1", vacuum_solid_LV, false, 0, fCheckOverlaps);
ColdPlate1_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));

//OneKPlate
G4VSolid* OneKPlate_S = new G4Tubs( "OneKPlate", zeroRadius, OneKPlate_Radius, (OneKPlate_Thickness /2.0), startAngle, spanningAngle);
G4LogicalVolume *OneKPlate_LV = new G4LogicalVolume( OneKPlate_S, Copper, "OneKPlate" );
OneKPlate_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0, 0, MixingPlate_z+ColdPlate2_z+ColdPlate1_z+OneKPlate_z),OneKPlate_LV, "OneKPlate", vacuum_solid_LV, false, 0, fCheckOverlaps);
OneKPlate_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));

G4double z[4]=  {-OneKShield_Width, 0.0, 0.0, OneKShield_Height};
G4double ri[4]= {0.0, 0.0, OneKShield_Radius ,  OneKShield_Radius };
G4double ro[4]= {OneKShield_Radius+OneKShield_Width, OneKShield_Radius+OneKShield_Width, OneKShield_Radius+OneKShield_Width, OneKShield_Radius+OneKShield_Width};
G4ThreeVector position_OneKShield = G4ThreeVector(0, 0, MixingPlate_z+ColdPlate2_z+ColdPlate1_z+OneKPlate_z) + G4ThreeVector(0, 0, -OneKPlate_Thickness/2-OneKShield_Height);

G4Polycone* OneKShield_S = new G4Polycone("OneKShield", startAngle, spanningAngle, 4, z, ri, ro);
G4LogicalVolume*  OneKShield_LV= new G4LogicalVolume(OneKShield_S, Copper, "OneKShield");
//G4LogicalVolume*  OneKShield_LV= new G4LogicalVolume(OneKShield_S, Aluminum, "OneKShield");
OneKShield_PV = new G4PVPlacement(NO_ROT, position_OneKShield , OneKShield_LV, "OneKShield", vacuum_solid_LV, false, 0, fCheckOverlaps);
OneKShield_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));

//FourKPlate
G4VSolid* FourKPlate_S = new G4Tubs( "FourKPlate", zeroRadius, FourKPlate_Radius, (FourKPlate_Thickness/2.0), startAngle, spanningAngle);
G4LogicalVolume *FourKPlate_LV = new G4LogicalVolume( FourKPlate_S, Copper, "FourKPlate" );
FourKPlate_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0, 0, MixingPlate_z+ColdPlate2_z+ColdPlate1_z+OneKPlate_z+FourKPlate_z),FourKPlate_LV, "FourKPlate", vacuum_solid_LV, false, 0, fCheckOverlaps);
FourKPlate_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));


G4double z1[4]=  {-FourKShield_Width, 0.0, 0.0, FourKShield_Height};
G4double ri1[4]= {0.0, 0.0, FourKShield_Radius ,  FourKShield_Radius };
G4double ro1[4]= {FourKShield_Radius+FourKShield_Width, FourKShield_Radius+FourKShield_Width, FourKShield_Radius+FourKShield_Width, FourKShield_Radius+FourKShield_Width};
G4ThreeVector position_FourKShield = G4ThreeVector(0, 0, MixingPlate_z+ColdPlate2_z+ColdPlate1_z+OneKPlate_z+FourKPlate_z) + G4ThreeVector(0, 0, -FourKPlate_Thickness/2-FourKShield_Height);

G4Polycone* FourKShield_S = new G4Polycone("FourKShield", startAngle, spanningAngle, 4, z1, ri1, ro1);
G4LogicalVolume*  FourKShield_LV= new G4LogicalVolume(FourKShield_S, Copper, "FourKShield");
//G4LogicalVolume*  FourKShield_LV= new G4LogicalVolume(FourKShield_S, Aluminum, "FourKShield");
FourKShield_PV = new G4PVPlacement(NO_ROT, position_FourKShield , FourKShield_LV, "FourKShield", vacuum_solid_LV, false, 0, fCheckOverlaps);
FourKShield_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));


//SeventyKPlate
G4VSolid* SeventyKPlate_S = new G4Tubs( "SeventyKPlate", zeroRadius, SeventyKPlate_Radius, (SeventyKPlate_Thickness /2.0), startAngle, spanningAngle);
G4LogicalVolume *SeventyKPlate_LV = new G4LogicalVolume( SeventyKPlate_S, Copper, "SeventyKPlate" );
SeventyKPlate_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0, 0, MixingPlate_z+ColdPlate2_z+ColdPlate1_z+OneKPlate_z+FourKPlate_z+SeventyKPlate_z),SeventyKPlate_LV, "SeventyKPlate", vacuum_solid_LV, false, 0, fCheckOverlaps);
SeventyKPlate_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));


G4double z2[4]=  {-SeventyKShield_Width, 0.0, 0.0, SeventyKShield_Height};
G4double ri2[4]= {0.0, 0.0, SeventyKShield_Radius ,  SeventyKShield_Radius };
G4double ro2[4]= {SeventyKShield_Radius+SeventyKShield_Width, SeventyKShield_Radius+SeventyKShield_Width, SeventyKShield_Radius+SeventyKShield_Width, SeventyKShield_Radius+SeventyKShield_Width};
G4ThreeVector position_SeventyKShield = G4ThreeVector(0, 0, MixingPlate_z+ColdPlate2_z+ColdPlate1_z+OneKPlate_z+FourKPlate_z+SeventyKPlate_z) + G4ThreeVector(0, 0, -SeventyKPlate_Thickness/2-SeventyKShield_Height);

G4Polycone* SeventyKShield_S = new G4Polycone("SeventyKShield", startAngle, spanningAngle, 4, z2, ri2, ro2);
G4LogicalVolume*  SeventyKShield_LV= new G4LogicalVolume(SeventyKShield_S, Copper, "SeventyKShield");
//G4LogicalVolume*  SeventyKShield_LV= new G4LogicalVolume(SeventyKShield_S, Aluminum, "SeventyKShield");
SeventyKShield_PV = new G4PVPlacement(NO_ROT, position_SeventyKShield , SeventyKShield_LV, "SeventyKShield", vacuum_solid_LV, false, 0, fCheckOverlaps);
SeventyKShield_LV->SetVisAttributes(G4VisAttributes(G4Colour(1.,0.5,0.)));;


//OVC
G4double z3[4]=  {-OVCShield_Width, 0.0, 0.0, OVCShield_Height};
G4double ri3[4]= {0.0, 0.0, OVCShield_Radius ,  OVCShield_Radius };
G4double ro3[4]= {OVCShield_Radius+OVCShield_Width, OVCShield_Radius+OVCShield_Width, OVCShield_Radius+OVCShield_Width, OVCShield_Radius+OVCShield_Width};
G4ThreeVector position_OVCShield = G4ThreeVector(0, 0, MixingPlate_z+ColdPlate2_z+ColdPlate1_z+OneKPlate_z+FourKPlate_z+SeventyKPlate_z+OVC_z) + G4ThreeVector(0, 0,-OVCShield_Height);

G4Polycone* OVCShield_S = new G4Polycone("OVCShield", startAngle, spanningAngle, 4, z3, ri3, ro3);
G4LogicalVolume*  OVCShield_LV= new G4LogicalVolume(OVCShield_S, Aluminum, "OVCShield");
OVCShield_PV = new G4PVPlacement(NO_ROT, position_OVCShield , OVCShield_LV, "OVCShield", vacuum_solid_LV, false, 0, fCheckOverlaps);
OVCShield_LV->SetVisAttributes(G4VisAttributes(G4Colour(G4Colour::Cyan())));
//OVCShield_LV->SetVisAttributes(G4VisAttributes::Invisible);



//Insulation but this is actually a surface to see the gammas coming out of the concrete walls and floor
//G4VSolid* Insulation_S = new G4Box("Insulation", Water_cylindercal_can_radius/2.0, delta/2.0, (Water_cylindercal_can_height+ConcreteSupport_height)/2.);
G4VSolid* Insulation_S = new G4Box("Insulation",lab68_wall_x/1.5, lab68_wall_y/1.5, lab68_wall_z);

//G4VSolid* test_2_S = new G4Box("test_2_solid", test_x/2.0, test_y/2.0 , test_z/2.0);
//G4VSolid* test_2_S = new G4Box("test_2_solid", test_x/2.0, test_y/2.0 , test_z/2.0);
//G4VSolid* test_hole_2_S = new G4Box("test_hole_2_solid", (test_x-2*delta)/2.0, (test_y-2*delta)/2.0, (test_z-2*delta)/2.0);
//G4SubtractionSolid* Insulation_S= new G4SubtractionSolid("Insulation", test_2_S, test_hole_2_S, NO_ROT, G4ThreeVector(0.,0., 0.));
G4LogicalVolume* Insulation_LV = new G4LogicalVolume(Insulation_S, Vacuum, "Insulation");
//Insulation_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(0, 0, 36.7*cm + 9.6*cm + 9.4*cm) , Insulation_LV, "Insulation", vacuum_solid_LV, false, 0, fCheckOverlaps);
//Insulation_PV = new G4PVPlacement(NO_ROT, G4ThreeVector{0.5*m,-1.2*m,0.75*m}, Insulation_LV, "Insulation", vacuum_solid_LV, false, 0, fCheckOverlaps);
////Insulation_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(0., fFilterCellSpacing-colimator_length-3*delta/2.0, (Water_cylindercal_can_height-ConcreteSupport_height)/2 - DT_Ti_T_location - Insulation_Thickness), Insulation_LV, "Insulation", vacuum_solid_LV, false, 0, fCheckOverlaps);
Insulation_LV->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));







//
// Visualization attributes
//
  vacuum_solid_LV->SetVisAttributes(G4VisAttributes::Invisible);
  //inner_shield_LV->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* test_vis = new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.5));  //used in verification of the geometry
  //G4VisAttributes* iron_vis = new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.5));
  G4VisAttributes* aluminum_vis = new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.5));


  //G4VisAttributes* lead_vis = new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.5));
  G4VisAttributes* helium_vis = new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.5));
  G4VisAttributes* lead_vis = new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.5));
  G4VisAttributes* silicon_vis = new G4VisAttributes(G4Colour(0.0,1.0,0.0,0.5));
  G4VisAttributes* ej_254_5pc_vis = new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.5));
  G4VisAttributes* LithiumIodide_vis = new G4VisAttributes(G4Colour(0.6,0.5,1.0,0.5));
  G4VisAttributes* air_vis = new G4VisAttributes(G4Colour(1.0,0.5,0.5,0.5));
  G4VisAttributes* EJ230_vis = new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.5));
  G4VisAttributes* EJ4265HD_vis = new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.5));

  //
  // Always return the physical World
  //
  return vacuum_solid_PV;
}
