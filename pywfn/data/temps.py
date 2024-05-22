"""
记录一些模板文件
"""

gjf=\
"""%chk=<CHK>.chk
# opt freq b3lyp/6-31g(d) gfinput pop=full iop(3/33=1)

Title Card Required

<CHARGE> <MULTI>
<COORD>

"""
xyz=\
"""<NATM>
<TITLE>
<COORD>
"""

si=\
"""<FILENAME>

<ENERGY>

Cartesian coordinate
<COORD>

Variable frequencies
<FREQ>
"""

pes=\
"""<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
<CDXML
 CreationProgram="ChemDraw 20.0.0.41"
 Name="none.cdxml"
 BoundingBox="0 0 0 0"
 WindowPosition="0 0"
 WindowSize="-2147483648 -2147483648"
 WindowIsZoomed="yes"
 FractionalWidths="yes"
 InterpretChemically="yes"
 ShowAtomQuery="yes"
 ShowAtomStereo="no"
 ShowAtomEnhancedStereo="yes"
 ShowAtomNumber="no"
 ShowResidueID="no"
 ShowBondQuery="yes"
 ShowBondRxn="yes"
 ShowBondStereo="no"
 ShowTerminalCarbonLabels="no"
 ShowNonTerminalCarbonLabels="no"
 HideImplicitHydrogens="no"
 LabelFont="3"
 LabelSize="10"
 LabelFace="96"
 CaptionFont="3"
 CaptionSize="10"
 HashSpacing="2.50"
 MarginWidth="1.60"
 LineWidth="0.60"
 BoldWidth="2"
 BondLength="14.40"
 BondSpacing="18"
 ChainAngle="120"
 LabelJustification="Auto"
 CaptionJustification="Left"
 AminoAcidTermini="HOH"
 ShowSequenceTermini="yes"
 ShowSequenceBonds="yes"
 ShowSequenceUnlinkedBranches="no"
 ResidueWrapCount="40"
 ResidueBlockCount="10"
 ResidueZigZag="yes"
 NumberResidueBlocks="no"
 PrintMargins="36 36 36 36"
 MacPrintInfo="0003000001200120000000000B6608A0FF84FF880BE309180367052703FC0002000001200120000000000B6608A0000100000064000000010001010100000001270F000100010000000000000000000000000002001901900000000000600000000000000000000100000000000000000000000000000000"
 ChemPropName=""
 ChemPropFormula="Chemical Formula: "
 ChemPropExactMass="Exact Mass: "
 ChemPropMolWt="Molecular Weight: "
 ChemPropMOverZ="m/z: "
 ChemPropAnalysis="Elemental Analysis: "
 ChemPropBoilingPt="Boiling Point: "
 ChemPropMeltingPt="Melting Point: "
 ChemPropCritTemp="Critical Temp: "
 ChemPropCritPres="Critical Pres: "
 ChemPropCritVol="Critical Vol: "
 ChemPropGibbs="Gibbs Energy: "
 ChemPropLogP="Log P: "
 ChemPropMR="MR: "
 ChemPropHenry="Henry&apos;s Law: "
 ChemPropEForm="Heat of Form: "
 ChemProptPSA="tPSA: "
 ChemPropCLogP="CLogP: "
 ChemPropCMR="CMR: "
 ChemPropLogS="LogS: "
 ChemPropPKa="pKa: "
 ChemPropID=""
 ChemPropFragmentLabel=""
 color="0"
 bgcolor="1"
 RxnAutonumberStart="1"
 RxnAutonumberConditions="no"
 RxnAutonumberStyle="Roman"
 RxnAutonumberFormat="(#)"
>
<colortable>
<color r="1" g="1" b="1"/>
<color r="0" g="0" b="0"/>
</colortable>
<fonttable>
<font id="1" charset="iso-8859-1" name="Arial"/>
</fonttable>
<page
 id="2"
 BoundingBox="0 0 540 719.75"
 HeaderPosition="36"
 FooterPosition="36"
 PrintTrimMarks="yes"
 HeightPages="1"
 WidthPages="1"
>
<fragment
 id="3"
 BoundingBox="[BOUND]"
 Z="6"
>
[NBS]
</fragment>
[TXS]
</page>
</CDXML>"""

mol="""[TITLE]

Created by pywfn
[NATOM][NBOND]  0  0  0  0  0  0  0  0  0    0
[ATOMBLOCK]
[BONDBLOCK]
M  END
"""