# LoopMixArticle
Scripts for the reclassification of uncertain observation project.
This README file was generated on 2021-02-02 by Emy Guilbault

GENERAL INFORMATION

1. Title of Dataset: Myxophies species location and environmental variables

2. Author Information
	A. Principal Investigator Contact Information
		Name: Emy Guilbault
		Institution: The University of Newcastle
		Address: University Dr, Callaghan NSW 2308
		Email: guilbaultemy@gmail.com

	B. Associate or Co-investigator Contact Information
		Name: Ian Renner
		Institution: The University of Newcastle
		Address: University Dr, Callaghan NSW 2308
		Email: ian.renner@newcastle.edu.au

	C. Alternate Contact Information
		Name: Michael Mahony
		Institution: The University of Newcastle
		Address: University Dr, Callaghan NSW 2308
		Email: michael.mahony@newcastle.edu.au

3. Date of data collection (single date, range, approximate date) : The data were downloaded from ALA (https://www.ala.org.au/) and from 
		Oza, A.U., Lovett, K.E., Williams, S.E. & Moritz, C. (2012) Recent speciation and limited phylogeographic structure in mixophyes frogs from the Australian wet tropics. Molecular Phylogenetics and Evolution, 62, 407–413.

4. Geographic location of data collection: The data coordinates in the files are in Universal Transverse Mercator (UTM), Australia.

5. Information about funding sources that supported the collection of the data: No funding were required, we analyzed open data from an article and an online database.


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC0 for the processed data, Rscript and environmental information file created.

2. Links to publications that cite or use the data: NA

3. Links to other publicly accessible locations of the data: The data presented here are also available at https://github.com/EmyGlblt/LoopMixArticle.

4. Links/relationships to ancillary data sets:  NA

5. Was data derived from another source? yes
	A. If yes, list source(s): The data were downloaded from ALA (https://www.ala.org.au/) and from 
		Oza, A.U., Lovett, K.E., Williams, S.E. & Moritz, C. (2012) Recent speciation and limited phylogeographic structure in mixophyes frogs from the australian wet tropics. Molecular Phylogenetics and Evolution, 62, 407–413.
	For the environmental covariates, we used shapefiles to derive our data. These shapefiles are freely avaiable at BBCVL, UC Davis Biogeo group and Bureau of Meteorology.

6. Recommended citation for this dataset: Guilbault, Emy; Renner, Ian; Mahony, Michael; Beh, Eric (2021), Myxophies species location and environmental variables, Dryad, Dataset, https://doi.org/10.5061/dryad.vx0k6djqw


DATA & FILE OVERVIEW

1. File List: 
NewRealdata_Jan21.R: R script to analyse the data describe below.		
MyxophiesDATA.csv: Species data file containing the species label (Myxophies Carbinensis, schevilli or Coggeri or Unknown) as well as the coordinates of the observation.
BandGrid5k_220118.shp: Environmental covariates information calculated at a grid window of dimension 5kx5k.

2. Relationship between files, if important: The script allows to use both the species locations and the environmental information for the grid window to produce the analyses of our article.

3. Additional related data collected that was not included in the current data package: Simulation data and script are available in the github page mentioned before.

4. Are there multiple versions of the dataset? yes
	A. If yes, name of file(s) that was updated: 
		i. Why was the file updated? The file was updated to be cleaned and make it available for publication.
		ii. When was the file updated? 02/20/2021


METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 
Species data and environmental information were downloaded from open database as mentioned above. We use ArcGIS to create a grid window and measure the environmental covariate information at the location of the grid points.

2. Methods for processing the data: 
We use only the species label that were known from Ozra (2012) and the ALA (species identified after 2006). The other ALA data were considered as unknown and to be relabelled through our algorithm. We converted the specie coordinates to UTM.

3. Instrument- or software-specific information needed to interpret the data: 
To use the data, the software R is to be used with the following packages: 
library(spatstat)
library(lattice)
library(sp)
library(maptools)
library(raster)
library(geostatsp)
library(rgdal)
library(caret)
library(gridExtra)
library(grid)
library(latticeExtra)

4. Standards and calibration information, if appropriate: NA

5. Environmental/experimental conditions: NA

6. Describe any quality-assurance procedures performed on the data: The data were reviewed with Michael Mahony, specialist of this species.

7. People involved with sample collection, processing, analysis and/or submission: Emy Guilbault, Ian Renner, Michale Mahony and Eric Beh.


DATA-SPECIFIC INFORMATION FOR: NewRealdata_Jan21.R
1. Number of variables: none

2. Number of cases/rows: 542

3. Variable List: none

4. Missing data codes: none

5. Specialized formats or other abbreviations used: NA


DATA-SPECIFIC INFORMATION FOR: MyxophiesDATA.csv
1. Number of variables: 3

2. Number of cases/rows: 455

3. Variable List: Species, X, Y

4. Missing data codes: none

5. Specialized formats or other abbreviations used: Mcarb=Myxophies carbinensis, Mcog=Myxophies coggeri, Msch=Myxophies schevilli



DATA-SPECIFIC INFORMATION FOR: BandGrid5k_220118.shp

1. Number of variables: 95

2. Number of cases/rows: 1630

3. Variable List: FID_Gridba	FID_Grid_1	Id	X	Y	FID_NorthB	Name	FID_AUS_ad	ID_0	ISO	NAME_0	OBJECTID_1	ISO3	NAME_ENGLI	NAME_ISO	NAME_FAO	NAME_LOCAL	NAME_OBSOL	NAME_VARIA	NAME_NONLA	NAME_FRENC	NAME_SPANI	NAME_RUSSI	NAME_ARABI	NAME_CHINE	WASPARTOF	CONTAINS	SOVEREIGN	ISO2	WWW	FIPS	ISON	VALIDFR	VALIDTO	POP2000	SQKM	POPSQKM	UNREGION1	UNREGION2	DEVELOPING	CIS	Transition	OECD	WBREGION	WBINCOME	WBDEBT	WBOTHER	CEEAC	CEMAC	CEPLG	COMESA	EAC	ECOWAS	IGAD	IOC	MRU	SACU	UEMOA	UMA	PALOP	PARTA	CACM	EurAsEC	Agadir	SAARC	ASEAN	NAFTA	GCC	CSN	CARICOM	EU	CAN	ACP	Landlocked	AOSIS	SIDS	Islands	LDC	AUS_alt	clim05UTM1	clim11utm1	clim13utm1	CLIMOND_18	clim26utm1	clim27utm1	CLIM34UTM1	bio06	NEAR_DIST	NEAR_FC	NEAR_DIS_1	CLIM05_1k1	CLIM06_1k	CLIM11_1k	CLIM13_1k	CLIM18_1k
However the Rscript subset this document to make use of the only variable that matters: OBJECTID_1, AUS_alt, CLIM05_1k1, CLIM06_1k, CLIM11_1k, CLIM13_1k, CLIM18_1k, NEAR_DIST, NEAR_DIS_1, X, Y

4. Missing data codes: blank

5. Specialized formats or other abbreviations used: NEAR_DIST=nearest distance to a stream, NEAR_DIS_1=nearest distance to the road
