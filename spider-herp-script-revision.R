############################################################################
############################################################################
############################################################################
########## Relationship between predatory spiders and their 	############# 
########## prey frogs and lizards, at La Selva, Costa Rica		#############
########## 		    -- B Folt, March 2020 --		     			      #############
##########          Re-examined April 2021                    #############
############################################################################
############################################################################
############################################################################

### 11 April 2021
### NOTE: This script was written with R version 3.5 and this version must be used. 
###       More recent versions of R operate with different rules and fail to 
###       perform tasks. I switched from R v. 4.0.2 back to v. 3.5 by navigating 
###       Tools > Global Options > Browsing to select the folder with v. 3.5
###       (Macintosh HD > Library > Frameworks > R.frameworks > Versions > 3.5).
###       This caused the data handling languages (Section I) to operate as intended.
###       I also was unable to load RPresence using R v. 4.0.2, and was only able to
###       load and use the package when I had reverted back to R v. 3.5.


# Clear the workspace
rm(list=ls())

# Set the working directory
setwd("/Users/brian/Dropbox/Auburn Ph.D. Spider predators/Predator-prey-models-GitHub")

datum = read.csv("captures.csv", header = TRUE)	# Load the datafile

### Load package 'reshape' to access data manipulation code 

library(reshape)
library(tidyr)

###########################################################################
############# I) Manipulate frog, lizard, and spider survey ###############
############# data and arthropod/litter data for analysis   ############### 
###########################################################################

######### Use for-loops to:
######### A) subset to different species (CRABRA, NORHUM, OOPPUM, spiders) 
######### 	and generate detection histories for species in each cell
######### B) organize food abundance covariates for occupancy models

datum = droplevels(subset(datum, Location != "NA"))	

######## A -- Create a vector of the focal species

species = c("CRABRA", "OOPPUM", "NORHUM", "ctenid") 

################ Initiate for-loop for above tasks ################ 

for (s in species){

herp = droplevels(subset(datum, SpeciesCode == s))

#herp = droplevels(subset(herp, herp$SVL != 1))	
# This line eliminates small ctenids (size class = 1)
# and only focuses on large indivduals of size classes 2+3
# HOWEVER, I commented it out, b/c it doesn't change the results.

## Use for-loops to tabulate the number of individuals 
## detected in each grid cell of each plot 
## during each of the surveys in each month

DetHist = matrix(, nrow = 1, ncol = 5)
detections = matrix(NA, nrow = 21, ncol = 5)

# Initiate the DETECTION HISTORY for-loop

for (i in 1:14){					# For each plot, 1:14
	plot = subset(herp, TreeNo == i)	# Subset to a plot

for (j in 1:9){						# For each month, 1:9
	month = subset(plot, as.numeric(Sesh) == j)	# Subset to a month

for (k in 1:3){					# For each survey occasion
	survey = subset(month, Occ == k)		# Subset to a survey
					
(tab = as.data.frame(table(survey$Cell)))	# Tabulate captures by cell

detections[,k] = tab$Freq			# Fill in matrix w/ detections

if(max(as.numeric(month$Sesh)) == 2) 	# If it's March with only two surveys,
	{detections[,3] = NA} 		# Then third survey column are NA

}						# End for-loop for survey occasions

detections[,4] = j				# Fill in Month column of matrix
detections[,5] = i				# Fill in TreeNo column of matrix

DetHist = rbind(DetHist, detections)	
	# Binds detection history for surveys
	# in a plot and month to a running tally
	# of all surveys among all plots/months

}}	# END DETECTION HISTORY for-loops for month and plot

DetHist = DetHist[-c(1),]
head(DetHist)
colnames(DetHist)=c("S1","S2","S3","Month","TreeNo")
DetHist2 = as.data.frame(DetHist)


## Reshape the data to format it for 'unmarked'
DetHist2 = melt(DetHist2, id=c("Month","TreeNo"))
colnames(DetHist2) = c("Month","TreeNo","Survey","Detection")
head(DetHist2)

# Create a vector for the Months & TreeNo
Month = DetHist2$Month
TreeNo = DetHist2$TreeNo

# Create a vector describing dry season (Feb, March, May) and wet (June-Nov)
Season = Month
Season[Season==2] = 1
Season[Season==3] = 1
Season[Season==4] = 2
Season[Season==5] = 2
Season[Season==6] = 2
Season[Season==7] = 2
Season[Season==8] = 2
Season[Season==9] = 2

# Create a new vector, CellNo, which identifies each cell within tree plot
x = c(1,10,11,12,13,14,15,16,17,18,19,2,20,21,3,4,5,6,7,8,9)
CellNo = as.data.frame(rep(x,378))
colnames(CellNo) = c("CellNo")

# Change all the counts (which can be greater than 1) to be 1
Detection = DetHist2$Detection
Detection[ Detection > 0 ] = 1

### If I want to build single-season models, then: 
# Reformat the surveys to be numeric 1:3 in each month
Survey = as.numeric(DetHist2$Survey)

# Bind these manipulated files back together to generate file for 'unmarked'
DetHist3 = cbind(TreeNo,Month,Season,CellNo,Survey,Detection)
DetHist3 = DetHist3[order(DetHist3$TreeNo, DetHist3$Month,
		 DetHist3$CellNo, DetHist3$Survey),]
head(DetHist3,20)

# Spread the surveys within months to wide format
DetHist4 = spread(DetHist3, Survey, Detection, fill=NA)

# Apply a unique name to each reformatted datafile
assign(paste0(s,"detections"), DetHist4) 

} # End the species loop

head(CRABRAdetections,10)	# E.g.,
tail(CRABRAdetections,10)	# E.g.,


######## B -- Create vectors describing arthropod abundance

arthropods = read.csv("arthropod-and-litter-data-v2.csv", header=TRUE)
# In this datafile, I increased the number of arthropod observations.
# For each cell corner with a sample, I copied it for
# each adjacent or katy-cornered cell. 
# Only issue is, after binding this with the herp observations, 
# some herp detection histories in months will be repeatedly used
# for cells that had two or more adjacent litter samples

bugs = arthropods[,c("Tree.Number","TreeSpecies3","Month3","GridCells",
		"Dry.Litter.Mass..g.","Acari","Araneae","Coleoptera",
		"Isopoda","Formicidae","Orthoptera")]
colnames(bugs) = c("TreeNo","Treatment","Month","CellNo","LitterMass",
			"Acari","Araneae","Coleoptera","Isopoda",
			"Formicidae","Orthoptera")
head(bugs)


### Merge the detection histories with the observational covariates
CRABRAdetections = merge(CRABRAdetections, bugs,
	 all.y=TRUE, by=c("TreeNo","Month","CellNo"))

OOPPUMdetections  = merge(OOPPUMdetections, bugs,
   all.y=TRUE, by=c("TreeNo","Month","CellNo"))

NORHUMdetections = merge(NORHUMdetections, bugs,
	 all.y=TRUE, by=c("TreeNo","Month","CellNo"))

cteniddetections = merge(cteniddetections, bugs,
	 all.y=TRUE, by=c("TreeNo","Month","CellNo"))


# Bind each frog and lizard dataframe with the ctenid one
CTENIDvsCRABRA = rbind(cteniddetections,CRABRAdetections)
CTENIDvsOOPPUM = rbind(cteniddetections,OOPPUMdetections)
CTENIDvsNORHUM = rbind(cteniddetections,NORHUMdetections)
head(CTENIDvsNORHUM, 20)

# Analyze the dataframes in R using RPresence

### Tabulate the number of surveys, number of detections for
### frogs, lizards, and spiders

length(na.omit(subset(CRABRAdetections, LitterMass != "NA")[,5:7])[,1]) # Number of sites
length(na.omit(subset(CRABRAdetections, LitterMass != "NA")[,5:7])[,1])*3 # Number of surveys
sum(na.omit(subset(CRABRAdetections, LitterMass != "NA")[,5:7])) # Number of observations
na.omit(subset(CRABRAdetections, LitterMass != "NA")) # The data

length(na.omit(subset(OOPPUMdetections, LitterMass != "NA")[,5:7])[,1]) # Number of sites
length(na.omit(subset(OOPPUMdetections, LitterMass != "NA")[,5:7])[,1])*3 # Number of surveys
sum(na.omit(subset(OOPPUMdetections, LitterMass != "NA")[,5:7])) # Number of observations

length(na.omit(subset(NORHUMdetections, LitterMass != "NA")[,5:7])[,1]) # Number of sites
length(na.omit(subset(NORHUMdetections, LitterMass != "NA")[,5:7])[,1])*3 # Number of surveys
sum(na.omit(subset(NORHUMdetections, LitterMass != "NA")[,5:7])) # Number of observations

length(na.omit(subset(cteniddetections, LitterMass != "NA")[,5:7])[,1]) # Number of sites
length(na.omit(subset(cteniddetections, LitterMass != "NA")[,5:7])[,1])*3 # Number of surveys
sum(na.omit(subset(cteniddetections, LitterMass != "NA")[,5:7])) # Number of observations



###########################################################################
############# II) Use the R library 'RPresence' to perform	  #############
#############     single-season two-species occupancy 	  	  #############
#############	    models for each species in R			          #############
###########################################################################

# Load the package
install.packages('/Users/brian/Desktop/RPresence_2.12.24.tar.gz')
library(RPresence)	 


####### i) set up the dataframe for each species

pairnames = c("CTENIDvsCRABRA","CTENIDvsOOPPUM","CTENIDvsNORHUM")
pairs = list(CTENIDvsCRABRA,CTENIDvsOOPPUM,CTENIDvsNORHUM)

for (i in 1:3){

twospecies = pairs[[i]]
twospecies = na.omit(twospecies)
twospecies = twospecies[,c(5:7,4,9:15)]
	# Ctenids are the first half of rows, vertebrate are the second half
	# Cols 1-3 are spider/herps detections, 4-11 are covariates 

# Specify the number of sites and no. of surveys
NoSurveys=ncol(twospecies[,1:3])         
NoSites=nrow(twospecies)

# Specify detection history
DetHist=matrix(as.integer(unlist(twospecies[,c(1:3)])),nrow=NoSites) 

# Specify site (cov1) and observation (cov2) covariates
#cov1=c(log(twospecies[5:11]+1),twospecies[4])	# Log-transformed to help w/ convergence
cov1=log(twospecies[4:11]+1)				# Log-transformed to help w/ convergence
cov2=NULL							# No observation covs considered here

# Create PAO file
PAO = createPao(DetHist,unitcov=cov1,survcov=cov2)
assign(paste0("PAO",pairnames[i]), PAO)

}


####### ii) set up the models    
####### and run them    

########################################## 
####### A ) Craugastor bransfordii ####### 
#######     vs. Ctenus sp.         ####### 
##########################################

PAO = PAOCTENIDvsCRABRA

### For this section, let's take the original six models from section II, 
### and add covariates for litter and bugs. This will test whether 
### predators or food/litter is more important for influencing
### occupancy and detection of frogs and lizards.

### CRABRA is a dietary generalist, but selectively eats 
### Acari, Araneae, Coleoptera, and Isopoda at proportions greater
### than they are available. So, let's taylor the models for these taxa.

### Original six models w/o coviarates
m1 = list(psi~SP,p~SP)
m2 = list(psi~SP,p~SP+INT_o)
m3 = list(psi~SP,p~SP+INT_o+INT_d)

m4  = list(psi~SP+INT,p~SP)
m5 = list(psi~SP+INT,p~SP+INT_o)
m6 = list(psi~SP+INT,p~SP+INT_o+INT_d)

### 9 models with litter covariates BUT no predator effect
m7 = list(psi~SP+LitterMass,p~SP)
m8 = list(psi~SP+LitterMass,p~SP+INT_o)
m9 = list(psi~SP+LitterMass,p~SP+INT_o+INT_d)

m10 = list(psi~SP,p~SP+LitterMass)
m11 = list(psi~SP,p~SP+INT_o+LitterMass)
m12 = list(psi~SP,p~SP+INT_o+INT_d+LitterMass)

m13 = list(psi~SP+LitterMass,p~SP+LitterMass)
m14 = list(psi~SP+LitterMass,p~SP+INT_o+LitterMass)
m15 = list(psi~SP+LitterMass,p~SP+INT_o+INT_d+LitterMass)

# 9 models with litter AND predator effect
m16 = list(psi~SP+INT+LitterMass,p~SP)
m17 = list(psi~SP+INT+LitterMass,p~SP+INT_o)
m18 = list(psi~SP+INT+LitterMass,p~SP+INT_o+INT_d)

m19 = list(psi~SP+INT,p~SP+LitterMass)
m20 = list(psi~SP+INT,p~SP+INT_o+LitterMass)
m21 = list(psi~SP+INT,p~SP+INT_o+INT_d+LitterMass)

m22 = list(psi~SP+INT+LitterMass,p~SP+LitterMass)
m23 = list(psi~SP+INT+LitterMass,p~SP+INT_o+LitterMass)
m24 = list(psi~SP+INT+LitterMass,p~SP+INT_o+INT_d+LitterMass)

### 9 models with food covariates BUT no predator effect
m25 = list(psi~SP+I(Acari+Araneae+Coleoptera+Isopoda),p~SP)
m26 = list(psi~SP+I(Acari+Araneae+Coleoptera+Isopoda),p~SP+INT_o)
m27 = list(psi~SP+I(Acari+Araneae+Coleoptera+Isopoda),p~SP+INT_o+INT_d)

m28 = list(psi~SP,p~SP+I(Acari+Araneae+Coleoptera+Isopoda))
m29 = list(psi~SP,p~SP+INT_o+I(Acari+Araneae+Coleoptera+Isopoda))
m30 = list(psi~SP,p~SP+INT_o+INT_d+I(Acari+Araneae+Coleoptera+Isopoda))

m31 = list(psi~SP+I(Acari+Araneae+Coleoptera+Isopoda),p~SP+I(Acari+Araneae+Coleoptera+Isopoda))
m32 = list(psi~SP+I(Acari+Araneae+Coleoptera+Isopoda),p~SP+INT_o+I(Acari+Araneae+Coleoptera+Isopoda))
m33 = list(psi~SP+I(Acari+Araneae+Coleoptera+Isopoda),p~SP+INT_o+INT_d+I(Acari+Araneae+Coleoptera+Isopoda))

# 9 models with food AND predator effect
m34 = list(psi~SP+INT+I(Acari+Araneae+Coleoptera+Isopoda),p~SP)
m35 = list(psi~SP+INT+I(Acari+Araneae+Coleoptera+Isopoda),p~SP+INT_o)
m36 = list(psi~SP+INT+I(Acari+Araneae+Coleoptera+Isopoda),p~SP+INT_o+INT_d)

m37 = list(psi~SP+INT,p~SP+I(Acari+Araneae+Coleoptera+Isopoda))
m38 = list(psi~SP+INT,p~SP+INT_o+I(Acari+Araneae+Coleoptera+Isopoda))
m39 = list(psi~SP+INT,p~SP+INT_o+INT_d+I(Acari+Araneae+Coleoptera+Isopoda))

m40 = list(psi~SP+INT+I(Acari+Araneae+Coleoptera+Isopoda),p~SP+I(Acari+Araneae+Coleoptera+Isopoda))
m41 = list(psi~SP+INT+I(Acari+Araneae+Coleoptera+Isopoda),p~SP+INT_o+I(Acari+Araneae+Coleoptera+Isopoda))
m42 = list(psi~SP+INT+I(Acari+Araneae+Coleoptera+Isopoda),p~SP+INT_o+INT_d+I(Acari+Araneae+Coleoptera+Isopoda))

# 9 models with seasonal effects, no predator effect

m43 = list(psi~SP+Season,p~SP)
m44 = list(psi~SP+Season,p~SP+INT_o)
m45 = list(psi~SP+Season,p~SP+INT_o+INT_d)

m46 = list(psi~SP,p~SP+Season)
m47 = list(psi~SP,p~SP+INT_o+Season)
m48 = list(psi~SP,p~SP+INT_o+INT_d+Season)

m49 = list(psi~SP+Season,p~SP+Season)
m50 = list(psi~SP+Season,p~SP+INT_o+Season)
m51 = list(psi~SP+Season,p~SP+INT_o+INT_d+Season)

# 9 models with litter AND predator effect
m52 = list(psi~SP+INT+Season,p~SP)
m53 = list(psi~SP+INT+Season,p~SP+INT_o)
m54 = list(psi~SP+INT+Season,p~SP+INT_o+INT_d)

m55 = list(psi~SP+INT,p~SP+Season)
m56 = list(psi~SP+INT,p~SP+INT_o+Season)
m57 = list(psi~SP+INT,p~SP+INT_o+INT_d+Season)

m58 = list(psi~SP+INT+Season,p~SP+Season)
m59 = list(psi~SP+INT+Season,p~SP+INT_o+Season)
m60 = list(psi~SP+INT+Season,p~SP+INT_o+INT_d+Season)


# Combine all the models into a list
ModelList = list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,
		m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
		m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,
		m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
		m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,
		m51,m52,m53,m54,m55,m56,m57,m58,m59,m60)

N = length(ModelList) # Specify the number of models

## Run the 60 models

ptm = proc.time()	### START A TIMER!

for (i in 1:N) {
	model = occMod(model=ModelList[[i]], data=PAO, type="so.2sp.1")
	assign(paste0("m",i), model)
}

proc.time() - ptm ### END THE TIMER, 60 models takes ~3 min

ModelsCB = list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,
		m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
		m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,
		m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
		m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,
		m51,m52,m53,m54,m55,m56,m57,m58,m59,m60)

ModRes = createAicTable(ModelsCB)
(CBtab = summary(ModRes))

(CBtopmods = subset(CBtab, wgt > 0.05))

# I note that some models have warning messages about the number of significant digits 
# reached during ML convergence; "Parameter esimates converged to approximately #.## signifcant digits."
# I considered these errors and found a situation where Jim Hines discussed this topic 
# on the PRESENCE forum here: http://www.phidot.org/forum/viewtopic.php?f=11&t=3491
# He says that as long as the number of significant digits is >=2, parameter estimates are accurate 
# to 3-4 decimal places. Since none of the top models with all the model weight have thrown the error 
# beneath this threshold, then my parameter estimates should be okay, and I am moving on.


# Copy and paste top-model set into Excel Table 2
clip = pipe("pbcopy","w")
write.table(CBtopmods, file=clip, sep="\t", row.names=FALSE)
close(clip)
# TABLE 2


###################################### 
####### B) Oophaga pumilio     ####### 
#######     vs. Ctenus sp.     ####### 
######################################

PAO = PAOCTENIDvsOOPPUM

### For this section, let's take the original eight models from section II, 
### and add covariates for litter and bugs. This will test whether 
### predators or food/litter is more important for influencing
### occupancy and detection of frogs and lizards.

### OOPPUM is a dietary specialist, which selectively eats 
### Acari and Formicidae. So, let's taylor the models for these taxa.

### Original six models w/o coviarates
m1 = list(psi~SP,p~SP)
m2 = list(psi~SP,p~SP+INT_o)
m3 = list(psi~SP,p~SP+INT_o+INT_d)

m4  = list(psi~SP+INT,p~SP)
m5 = list(psi~SP+INT,p~SP+INT_o)
m6 = list(psi~SP+INT,p~SP+INT_o+INT_d)

### 9 models with litter covariates BUT no predator effect
m7 = list(psi~SP+LitterMass,p~SP)
m8 = list(psi~SP+LitterMass,p~SP+INT_o)
m9 = list(psi~SP+LitterMass,p~SP+INT_o+INT_d)

m10 = list(psi~SP,p~SP+LitterMass)
m11 = list(psi~SP,p~SP+INT_o+LitterMass)
m12 = list(psi~SP,p~SP+INT_o+INT_d+LitterMass)

m13 = list(psi~SP+LitterMass,p~SP+LitterMass)
m14 = list(psi~SP+LitterMass,p~SP+INT_o+LitterMass)
m15 = list(psi~SP+LitterMass,p~SP+INT_o+INT_d+LitterMass)

# 9 models with litter AND predator effect
m16 = list(psi~SP+INT+LitterMass,p~SP)
m17 = list(psi~SP+INT+LitterMass,p~SP+INT_o)
m18 = list(psi~SP+INT+LitterMass,p~SP+INT_o+INT_d)

m19 = list(psi~SP+INT,p~SP+LitterMass)
m20 = list(psi~SP+INT,p~SP+INT_o+LitterMass)
m21 = list(psi~SP+INT,p~SP+INT_o+INT_d+LitterMass)

m22 = list(psi~SP+INT+LitterMass,p~SP+LitterMass)
m23 = list(psi~SP+INT+LitterMass,p~SP+INT_o+LitterMass)
m24 = list(psi~SP+INT+LitterMass,p~SP+INT_o+INT_d+LitterMass)

### 9 models with food covariates BUT no predator effect
m25 = list(psi~SP+I(Acari+Formicidae),p~SP)
m26 = list(psi~SP+I(Acari+Formicidae),p~SP+INT_o)
m27 = list(psi~SP+I(Acari+Formicidae),p~SP+INT_o+INT_d)

m28 = list(psi~SP,p~SP+I(Acari+Formicidae))
m29 = list(psi~SP,p~SP+INT_o+I(Acari+Formicidae))
m30 = list(psi~SP,p~SP+INT_o+INT_d+I(Acari+Formicidae))

m31 = list(psi~SP+I(Acari+Formicidae),p~SP+I(Acari+Formicidae))
m32 = list(psi~SP+I(Acari+Formicidae),p~SP+INT_o+I(Acari+Formicidae))
m33 = list(psi~SP+I(Acari+Formicidae),p~SP+INT_o+INT_d+I(Acari+Formicidae))

# 9 models with food AND predator effect
m34 = list(psi~SP+INT+I(Acari+Formicidae),p~SP)
m35 = list(psi~SP+INT+I(Acari+Formicidae),p~SP+INT_o)
m36 = list(psi~SP+INT+I(Acari+Formicidae),p~SP+INT_o+INT_d)

m37 = list(psi~SP+INT,p~SP+I(Acari+Formicidae))
m38 = list(psi~SP+INT,p~SP+INT_o+I(Acari+Formicidae))
m39 = list(psi~SP+INT,p~SP+INT_o+INT_d+I(Acari+Formicidae))

m40 = list(psi~SP+INT+I(Acari+Formicidae),p~SP+I(Acari+Formicidae))
m41 = list(psi~SP+INT+I(Acari+Formicidae),p~SP+INT_o+I(Acari+Formicidae))
m42 = list(psi~SP+INT+I(Acari+Formicidae),p~SP+INT_o+INT_d+I(Acari+Formicidae))

# 9 models with seasonal effects, no predator effect

m43 = list(psi~SP+Season,p~SP)
m44 = list(psi~SP+Season,p~SP+INT_o)
m45 = list(psi~SP+Season,p~SP+INT_o+INT_d)

m46 = list(psi~SP,p~SP+Season)
m47 = list(psi~SP,p~SP+INT_o+Season)
m48 = list(psi~SP,p~SP+INT_o+INT_d+Season)

m49 = list(psi~SP+Season,p~SP+Season)
m50 = list(psi~SP+Season,p~SP+INT_o+Season)
m51 = list(psi~SP+Season,p~SP+INT_o+INT_d+Season)

# 9 models with litter AND predator effect
m52 = list(psi~SP+INT+Season,p~SP)
m53 = list(psi~SP+INT+Season,p~SP+INT_o)
m54 = list(psi~SP+INT+Season,p~SP+INT_o+INT_d)

m55 = list(psi~SP+INT,p~SP+Season)
m56 = list(psi~SP+INT,p~SP+INT_o+Season)
m57 = list(psi~SP+INT,p~SP+INT_o+INT_d+Season)

m58 = list(psi~SP+INT+Season,p~SP+Season)
m59 = list(psi~SP+INT+Season,p~SP+INT_o+Season)
m60 = list(psi~SP+INT+Season,p~SP+INT_o+INT_d+Season)


# Combine all the models into a list
ModelList = list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,
		m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
		m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,
		m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
		m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,
		m51,m52,m53,m54,m55,m56,m57,m58,m59,m60)

N = length(ModelList) # Specify the number of models

## Run the 60 models

ptm = proc.time()	### START A TIMER!

for (i in 1:N) {
	model = occMod(model=ModelList[[i]], data=PAO, type="so.2sp.1")
	assign(paste0("m",i), model)
}

proc.time() - ptm ### END THE TIMER, 60 models takes ~3 min

ModelsOP = list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,
		m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
		m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,
		m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
		m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,
		m51,m52,m53,m54,m55,m56,m57,m58,m59,m60)

ModRes = createAicTable(ModelsOP)
(OPtab = summary(ModRes))

(OPtopmods = subset(OPtab, wgt > 0.05))

# Copy and paste top-model set into Excel Table 2
clip = pipe("pbcopy","w")
write.table(OPtopmods, file=clip, sep="\t", row.names=FALSE)
close(clip)
# TABLE 2



##################################### 
####### C) Norops humilis     ####### 
#######    vs. Ctenus sp.     ####### 
#####################################

PAO = PAOCTENIDvsNORHUM

### For this section, let's take the original eight models from section II, 
### and add covariates for litter and bugs. This will test whether 
### predators or food/litter is more important for influencing
### occupancy and detection of frogs and lizards.

### Original six models w/o coviarates
m1 = list(psi~SP,p~SP)
m2 = list(psi~SP,p~SP+INT_o)
m3 = list(psi~SP,p~SP+INT_o+INT_d)

m4  = list(psi~SP+INT,p~SP)
m5 = list(psi~SP+INT,p~SP+INT_o)
m6 = list(psi~SP+INT,p~SP+INT_o+INT_d)

### 9 models with litter covariates BUT no predator effect
m7 = list(psi~SP+LitterMass,p~SP)
m8 = list(psi~SP+LitterMass,p~SP+INT_o)
m9 = list(psi~SP+LitterMass,p~SP+INT_o+INT_d)

m10 = list(psi~SP,p~SP+LitterMass)
m11 = list(psi~SP,p~SP+INT_o+LitterMass)
m12 = list(psi~SP,p~SP+INT_o+INT_d+LitterMass)

m13 = list(psi~SP+LitterMass,p~SP+LitterMass)
m14 = list(psi~SP+LitterMass,p~SP+INT_o+LitterMass)
m15 = list(psi~SP+LitterMass,p~SP+INT_o+INT_d+LitterMass)

# 9 models with litter AND predator effect
m16 = list(psi~SP+INT+LitterMass,p~SP)
m17 = list(psi~SP+INT+LitterMass,p~SP+INT_o)
m18 = list(psi~SP+INT+LitterMass,p~SP+INT_o+INT_d)

m19 = list(psi~SP+INT,p~SP+LitterMass)
m20 = list(psi~SP+INT,p~SP+INT_o+LitterMass)
m21 = list(psi~SP+INT,p~SP+INT_o+INT_d+LitterMass)

m22 = list(psi~SP+INT+LitterMass,p~SP+LitterMass)
m23 = list(psi~SP+INT+LitterMass,p~SP+INT_o+LitterMass)
m24 = list(psi~SP+INT+LitterMass,p~SP+INT_o+INT_d+LitterMass)

### 9 models with food covariates BUT no predator effect
m25 = list(psi~SP+I(Araneae+Isopoda),p~SP)
m26 = list(psi~SP+I(Araneae+Isopoda),p~SP+INT_o)
m27 = list(psi~SP+I(Araneae+Isopoda),p~SP+INT_o+INT_d)

m28 = list(psi~SP,p~SP+I(Araneae+Isopoda))
m29 = list(psi~SP,p~SP+INT_o+I(Araneae+Isopoda))
m30 = list(psi~SP,p~SP+INT_o+INT_d+I(Araneae+Isopoda))

m31 = list(psi~SP+I(Araneae+Isopoda),p~SP+I(Araneae+Isopoda))
m32 = list(psi~SP+I(Araneae+Isopoda),p~SP+INT_o+I(Araneae+Isopoda))
m33 = list(psi~SP+I(Araneae+Isopoda),p~SP+INT_o+INT_d+I(Araneae+Isopoda))

# 9 models with food AND predator effect
m34 = list(psi~SP+INT+I(Araneae+Isopoda),p~SP)
m35 = list(psi~SP+INT+I(Araneae+Isopoda),p~SP+INT_o)
m36 = list(psi~SP+INT+I(Araneae+Isopoda),p~SP+INT_o+INT_d)

m37 = list(psi~SP+INT,p~SP+I(Araneae+Isopoda))
m38 = list(psi~SP+INT,p~SP+INT_o+I(Araneae+Isopoda))
m39 = list(psi~SP+INT,p~SP+INT_o+INT_d+I(Araneae+Isopoda))

m40 = list(psi~SP+INT+I(Araneae+Isopoda),p~SP+I(Araneae+Isopoda))
m41 = list(psi~SP+INT+I(Araneae+Isopoda),p~SP+INT_o+I(Araneae+Isopoda))
m42 = list(psi~SP+INT+I(Araneae+Isopoda),p~SP+INT_o+INT_d+I(Araneae+Isopoda))

# 9 models with seasonal effects, no predator effect

m43 = list(psi~SP+Season,p~SP)
m44 = list(psi~SP+Season,p~SP+INT_o)
m45 = list(psi~SP+Season,p~SP+INT_o+INT_d)

m46 = list(psi~SP,p~SP+Season)
m47 = list(psi~SP,p~SP+INT_o+Season)
m48 = list(psi~SP,p~SP+INT_o+INT_d+Season)

m49 = list(psi~SP+Season,p~SP+Season)
m50 = list(psi~SP+Season,p~SP+INT_o+Season)
m51 = list(psi~SP+Season,p~SP+INT_o+INT_d+Season)

# 9 models with litter AND predator effect
m52 = list(psi~SP+INT+Season,p~SP)
m53 = list(psi~SP+INT+Season,p~SP+INT_o)
m54 = list(psi~SP+INT+Season,p~SP+INT_o+INT_d)

m55 = list(psi~SP+INT,p~SP+Season)
m56 = list(psi~SP+INT,p~SP+INT_o+Season)
m57 = list(psi~SP+INT,p~SP+INT_o+INT_d+Season)

m58 = list(psi~SP+INT+Season,p~SP+Season)
m59 = list(psi~SP+INT+Season,p~SP+INT_o+Season)
m60 = list(psi~SP+INT+Season,p~SP+INT_o+INT_d+Season)

# Combine all the models into a list
ModelList = list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,
		m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
		m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,
		m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
		m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,
		m51,m52,m53,m54,m55,m56,m57,m58,m59,m60)

N = length(ModelList) # Specify the number of models

## Run the 60 models

ptm = proc.time()	### START A TIMER!

for (i in 1:N) {
	model = occMod(model=ModelList[[i]], data=PAO, type="so.2sp.1")
	assign(paste0("m",i), model)
}

proc.time() - ptm ### END THE TIMER, 60 models takes ~3 min

ModelsNH = list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,
		m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
		m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,
		m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
		m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,
		m51,m52,m53,m54,m55,m56,m57,m58,m59,m60)

ModRes = createAicTable(ModelsNH)
(NHtab = summary(ModRes))

(NHtopmods = subset(NHtab, wgt > 0.05))


# Copy and paste top-model set into Excel Table 2
clip = pipe("pbcopy","w")
write.table(NHtopmods, file=clip, sep="\t", row.names=FALSE)
close(clip)
# TABLE 2


## Now, need to average the parameters from across all the models 
## to get model-averaged parameter estimates. However, function 'mod.avg'
## doesn't work for two-species models. 

## So, I'll have to do this this process manually



###########################################################################
############# III) Calculate an average model for two-species #############
#############     occupancy and detection for CRABRA, OOPPUM, #############
#############	and NORHUM here MANUALLY 			  		  #############
###########################################################################

########################
####### A) CRABRA ######
########################

#### The model average code here RUNS, but i'm not sure it WORKS. Need to de-bug 
#### to make sure that the correct parameters are being stored through here. 

Models = ModelsCB

#### Occupancy

## Use for-loops to:
## model-average the coefficients and standard errors

## 1) Extract model name, number of parameters, AIC, and
##    parameter estimates and parameter SEs into a matrix

modlist = matrix(, nrow=N, ncol=15)

for (i in 1:N){

modlist[i,1] = Models[[i]][1]$modname
modlist[i,2] = Models[[i]][7]$npar
modlist[i,3] = Models[[i]][8]$aic
modlist[i,4] = Models[[i]][9]$beta$psi[1,1]
modlist[i,5] = Models[[i]][9]$beta$psi[2,1]

# Use if-else loops to fill in parameter values for each model
if(length(Models[[i]][9]$beta$psi[,1]) == 2) {		# 2 psi params
	modlist[i,6:9] = NA
} else {

if(length(Models[[i]][9]$beta$psi[,1]) == 3) {		# 3 psi params
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiBa"){
		modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # 3A; predator effect on psiA 
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiA.LitterMass_psiA"){
		modlist[i,7] = Models[[i]][9]$beta$psi[3,1] # 3B; litter effect on psiA
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiA.I(Acari_+_Araneae_+_Coleoptera_+_Isopoda)_psiA"){
		modlist[i,8] = Models[[i]][9]$beta$psi[3,1] # 3C; food effect on psiA
		} else {
		modlist[i,9] = Models[[i]][9]$beta$psi[3,1] # 3D; season effect on psiA
		}}}
} else {

if(length(Models[[i]][9]$beta$psi[,1]) == 4) {		# 4 psi params
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.LitterMass_psiA"){
	modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # Predator effect on psiA
	modlist[i,7] = Models[[i]][9]$beta$psi[4,1] # Litter effect on psiA
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.I(Acari_+_Araneae_+_Coleoptera_+_Isopoda)_psiA"){
	modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # Predator effect on psiA
	modlist[i,8] = Models[[i]][9]$beta$psi[4,1] # Food effect on psiA
		} else {
	modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # Predator effect on psiA
	modlist[i,9] = Models[[i]][9]$beta$psi[4,1] # Seasonal effect on psiA
    }
   }
  }
 }
}	# End parameter if-else loops

modlist[i,10] = Models[[i]][9]$beta$psi[1,2]
modlist[i,11] = Models[[i]][9]$beta$psi[2,2]



# Use if-else loops to fill in SE values for each model
if(length(Models[[i]][9]$beta$psi[,2]) == 2) {		# 2 psi params
	modlist[i,12:15] = NA
} else {

if(length(Models[[i]][9]$beta$psi[,1]) == 3) {		# 3 psi params
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiBa"){
		modlist[i,12] = Models[[i]][9]$beta$psi[3,2] # SE of predator effect
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiA.LitterMass_psiA"){
		modlist[i,13] = Models[[i]][9]$beta$psi[3,2] # SE of litter effect
			} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A4_psiA.I(Acari_+_Araneae_+_Coleoptera_+_Isopoda)_psiA"){
		modlist[i,14] = Models[[i]][9]$beta$psi[3,2] # SE of food effect
			} else {
		modlist[i,15] = Models[[i]][9]$beta$psi[3,2] # SE of seasonal effect
		}}}
} else {
	
if(length(Models[[i]][9]$beta$psi[,1]) == 4) {		# 4 psi params
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.LitterMass_psiA"){
	modlist[i,12] = Models[[i]][9]$beta$psi[3,2] # SE of predator effect
	modlist[i,13] = Models[[i]][9]$beta$psi[4,2] # SE of litter effect
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.I(Acari_+_Araneae_+_Coleoptera_+_Isopoda)_psiA"){
	modlist[i,12] = Models[[i]][9]$beta$psi[3,2] # SE of predator effect
	modlist[i,14] = Models[[i]][9]$beta$psi[4,2] # SE of food effect
		} else {
	modlist[i,12] = Models[[i]][9]$beta$psi[3,2] 
	modlist[i,15] = Models[[i]][9]$beta$psi[4,2] # SE of seasonal effect
    }
   }
  }
 }
}	# End SE if-else loops

}	# End model loop

modlist = as.data.frame(modlist)
head(modlist)

colnames(modlist) = c("ModelName","Pars","AIC",
			"Int","SP","INT2","Litter","Arthropods","Season",
			"seINT","seSP","seINT2","seLitter","seArthropods",
			"seSeason")

modlist[,c("Pars","AIC","Int","SP","INT2","Litter","Arthropods","Season",
	"seINT","seSP","seINT2","seLitter","seArthropods","seSeason")] = 
	lapply(modlist[,c("Pars","AIC","Int","SP","INT2","Litter",
	"Arthropods","Season","seINT","seSP","seINT2","seLitter",
	"seArthropods","seSeason")],
  	function(x) as.numeric(levels(x))[x])
head(modlist)

## 2) Calculate model weights  
##	for each model

modlist = modlist[order(modlist$AIC),]		# Sort by AIC

Sample = 872			# Sample size of sites
AICc = modlist$AIC+((2*modlist$Pars)*(modlist$Pars+1)/(Sample-modlist$Pars-1))
				# AICc; AIC corrected for small sample size
DeltaAICc = AICc - AICc[1]	# DeltaAICc
LL = exp(-0.5*DeltaAICc)	# Log-likelihood
Weight =LL/sum(LL)		# Model weight

ModRes = cbind(modlist[,1:3],AICc,DeltaAICc,LL,Weight,modlist[,4:15])


## 3) Calculate model-averaged coefficients estimates,
##	parameter weights, and unconditional standard error for each parameter!

AvgMod = matrix(NA,3,6, 
	dimnames=list(c("Coefficient","Parameter weight","Unconditional SE"),
	c("Intercept","SP","INT2","Litter","Arthropods","Season")))

subsets = list(c("Int","seINT","Weight"),c("SP","seSP","Weight"),
	c("INT2","seINT2","Weight"),c("Litter","seLitter","Weight"),
	c("Arthropods","seArthropods","Weight"),c("Season","seSeason","Weight"))
n = length(subsets)

for(i in 1:n){			# for each 'subset' from 1-6
	sub = ModRes[,subsets[[i]]]
	AvgMod[1,i] = sum(sub[,1]*sub$Weight,na.rm=TRUE)# Coefficient
	AvgMod[2,i] = sum(subset(sub,sub[,1] != "NA")[,3])# Param weight
	se = (sub[,2]^2+(sub[,1]-AvgMod[1,i])^2)^0.5
	se[is.na(se)] = AvgMod[1,i]
	bind = cbind(sub,se)
	SEsub = subset(bind, bind$se < summary(bind$se)[5])# Remove wild SEs
	AvgMod[3,i] = crossprod(SEsub$se, SEsub$Weight)}	# Unconditional SE
(AvgModCBpsi = AvgMod)


# Copy and paste the detection table into Excel for final manipulation
clip = pipe("pbcopy","w")
write.table(AvgModCBpsi, file=clip, sep="\t", row.names=FALSE)
close(clip)
# TABLE 3




##### DETECTION ESTIMATES 

## Use for-loops to model-average coefficients for DETECTION

## 1) Extract model name, number of parameters, AIC, and
##    parameter estimates and parameter SEs into a matrix

modlist = matrix(, nrow=N, ncol=17)

for (i in 1:N){		# Model loop

modlist[i,1] = Models[[i]][1]$modname
modlist[i,2] = Models[[i]][7]$npar
modlist[i,3] = Models[[i]][8]$aic
modlist[i,4] = Models[[i]][9]$beta$p[1,1]
modlist[i,5] = Models[[i]][9]$beta$p[2,1]

# Use if-else loops to fill in parameter values for each model
if(length(Models[[i]][9]$beta$p[,1]) == 2) {		# 2 p params
	modlist[i,6:10] = NA
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 3) {		# 3 p params
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_rA[1]"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] # rBA effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].LitterMass_pA"){
		modlist[i,8] = Models[[i]][9]$beta$p[3,1] # Litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].I(Acari_+_Araneae_+_Coleoptera_+_Isopoda)_pA"){
		modlist[i,9] = Models[[i]][9]$beta$p[3,1] # Food effect on detection
	} else {		
		modlist[i,10] = Models[[i]][9]$beta$p[3,1] # Seasonal effect on detection
	}}}
	
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 4) {		# 4 p params
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_rBa[1]"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] # rBA effect on detection
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] # rBa effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].LitterMass_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,8] = Models[[i]][9]$beta$p[4,1] # Litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].I(Acari_+_Araneae_+_Coleoptera_+_Isopoda)_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,9] = Models[[i]][9]$beta$p[4,1] # Food effect on detection
	} else {	
		modlist[i,6] = Models[[i]][9]$beta$p[3,1]  
		modlist[i,10] = Models[[i]][9]$beta$p[4,1] # Seasonal effect on detection
	}}}
	
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 5) {		# 5 p params
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].LitterMass_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] # rBA effect on detection
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] # rBa effect on detection
		modlist[i,8] = Models[[i]][9]$beta$p[5,1] # Litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].I(Acari_+_Araneae_+_Coleoptera_+_Isopoda)_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] 
		modlist[i,9] = Models[[i]][9]$beta$p[5,1] # Food effect on detection
	} else {	
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] 
		modlist[i,10] = Models[[i]][9]$beta$p[5,1] # Seasonal effect on detection
	 }
    }
   }
  }
 }
}	# End parameter if-else loops


modlist[i,11] = Models[[i]][9]$beta$p[1,2]
modlist[i,12] = Models[[i]][9]$beta$p[2,2]

# Use if-else loops to fill in SE values for each model
if(length(Models[[i]][9]$beta$p[,1]) == 2) {		# 2 p params
	modlist[i,13:16] = NA
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 3) {		# 3 p params
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_rA[1]"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] # SE rBA effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].LitterMass_pA"){
		modlist[i,15] = Models[[i]][9]$beta$p[3,2] # SE litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].I(Acari_+_Araneae_+_Coleoptera_+_Isopoda)_pA"){
		modlist[i,16] = Models[[i]][9]$beta$p[3,2] # SE food effect on detection
	} else {	
		modlist[i,17] = Models[[i]][9]$beta$p[3,2] # SE seasonal effect on detection
	}}}
	
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 4) {		# 4 p params
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_rBa[1]"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] # SE of rBA effect on detection
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] # SE of rBa effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].LitterMass_pA"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] # 
		modlist[i,15] = Models[[i]][9]$beta$p[4,2] # SE of litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].I(Acari_+_Araneae_+_Coleoptera_+_Isopoda)_pA"){
		modlist[i,14] = Models[[i]][9]$beta$p[3,2] # 
		modlist[i,16] = Models[[i]][9]$beta$p[4,2] # SE of food effect on detection
	} else {
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,17] = Models[[i]][9]$beta$p[4,2] # SE of  - Arthropods
	}}}
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 5) {		# 5 p params
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].LitterMass_pA"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] 
		modlist[i,15] = Models[[i]][9]$beta$p[5,2] # SE of litter effect
	} else {
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].I(Acari_+_Araneae_+_Coleoptera_+_Isopoda)_pA"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] 
		modlist[i,16] = Models[[i]][9]$beta$p[5,2] # SE of food effect				
	} else {
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] 
		modlist[i,17] = Models[[i]][9]$beta$p[5,2] # SE of season effect
	 }
    }
   }
  }
 }
}	# End SE if-else loops

}	# End model loop

modlist = as.data.frame(modlist)
head(modlist)

colnames(modlist) = c("ModelName","Pars","AIC",
			"Int","SP2","INTo","INTd","Litter","Arthropods","Season",
			"seINT","seSP2","seINTo","seINTd","seLitter","seArthropods","seSeason")

modlist[,c("Pars","AIC","Int","SP2","INTo","INTd","Litter","Arthropods","Season",
	"seINT","seSP2","seINTo","seINTd","seLitter","seArthropods","seSeason")] = 
	lapply(modlist[,c("Pars","AIC","Int","SP2","INTo","INTd","Litter",
	"Arthropods","Season","seINT","seSP2","seINTo","seINTd","seLitter",
	"seArthropods","seSeason")],
  	function(x) as.numeric(levels(x))[x])
head(modlist)

## 2) Calculate model weights  
##	for each model

modlist = modlist[order(modlist$AIC),]		# Sort by AIC

Sample = 872			# Sample size of sites
AICc = modlist$AIC+((2*modlist$Pars)*(modlist$Pars+1)/(Sample-modlist$Pars-1))
				# AICc; AIC corrected for small sample size
DeltaAICc = AICc - AICc[1]	# DeltaAICc
LL = exp(-0.5*DeltaAICc)	# Log-likelihood
Weight = LL/sum(LL)		# Model weight

ModRes = cbind(modlist[,1:3],AICc,DeltaAICc,LL,Weight,modlist[,4:17])
head(ModRes)

## 3) Calculate model-averaged coefficients estimates,
##	parameter weights, and unconditional standard error for each parameter!

AvgMod = matrix(NA,3,7, 
	dimnames=list(c("Coefficient","Parameter weight","Unconditional SE"),
	c("Intercept","SP2","INTo","INTd","Litter","Arthropods","Season")))

subsets = list(c("Int","seINT","Weight"),c("SP2","seSP2","Weight"),
	c("INTo","seINTo","Weight"),c("INTd","seINTd","Weight"),
	c("Litter","seLitter","Weight"),c("Arthropods","seArthropods","Weight"),
	c("Season","seSeason","Weight"))
n = length(subsets)

for(i in 1:n){
	sub = ModRes[,subsets[[i]]]
	AvgMod[1,i] = sum(sub[,1]*sub$Weight,na.rm=TRUE)	# Coefficient
	AvgMod[2,i] = sum(subset(sub,sub[,1] != "NA")[,3]) # Param weight
	se = (sub[,2]^2+(sub[,1]-AvgMod[1,i])^2)^0.5
	se[is.na(se)] = AvgMod[1,i]
	AvgMod[3,i] = crossprod(se,sub$Weight)}				# Unconditional SE
(AvgModCBp = AvgMod)

## Here I generated model-averaged DETECTION parameter estimates,
## parameter weights, and unconditional SEs.

## This information can be used to create:
## 1) A Table ('AvgMod') for the paper with the average model and param weights
## 2) Figures to show how detection varies conditionally, by covariates, or BOTH 

# Copy and paste the detection table into Excel for final manipulation
clip = pipe("pbcopy","w")
write.table(AvgModCBp, file=clip, sep="\t", row.names=FALSE)
close(clip)
# TABLE 4


########################
####### B) OOPPUM ######
######################## 

Models = ModelsOP

## Use for-loops to:
## model-average the coefficients and standard errors

#####
##### OCCUPANCY 
##### ESTIMATES
##### 

## 1) Extract model name, number of parameters, AIC, and
##    parameter estimates and parameter SEs into a matrix

modlist = matrix(, nrow=N, ncol=15)

for (i in 1:N){

modlist[i,1] = Models[[i]][1]$modname
modlist[i,2] = Models[[i]][7]$npar
modlist[i,3] = Models[[i]][8]$aic
modlist[i,4] = Models[[i]][9]$beta$psi[1,1]
modlist[i,5] = Models[[i]][9]$beta$psi[2,1]

# Use if-else loops to fill in parameter values for each model
if(length(Models[[i]][9]$beta$psi[,1]) == 2) {		# 2 psi params
	modlist[i,6:9] = NA
} else {

if(length(Models[[i]][9]$beta$psi[,1]) == 3) {		# 3 psi params
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiBa"){
		modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # 3A; predator effect on psiA 
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiA.LitterMass_psiA"){
		modlist[i,7] = Models[[i]][9]$beta$psi[3,1] # 3B; litter effect on psiA
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiA.I(Acari + Formicidae)_psiA"){
		modlist[i,8] = Models[[i]][9]$beta$psi[3,1] # 3C; food effect on psiA
		} else {
		modlist[i,9] = Models[[i]][9]$beta$psi[3,1] # 3D; season effect on psiA
		}}}
} else {

if(length(Models[[i]][9]$beta$psi[,1]) == 4) {		# 4 psi params
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.LitterMass_psiA"){
	modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # Predator effect on psiA
	modlist[i,7] = Models[[i]][9]$beta$psi[4,1] # Litter effect on psiA
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.I(Acari + Formicidae)_psiA"){
	modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # Predator effect on psiA
	modlist[i,8] = Models[[i]][9]$beta$psi[4,1] # Food effect on psiA
		} else {
	modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # Predator effect on psiA
	modlist[i,9] = Models[[i]][9]$beta$psi[4,1] # Seasonal effect on psiA
    }
   }
  }
 }
}	# End parameter if-else loops


# Insert SE for pA and pB
modlist[i,10] = Models[[i]][9]$beta$psi[1,2]
modlist[i,11] = Models[[i]][9]$beta$psi[2,2]

# Use if-else loops to fill in SE values for each model
if(length(Models[[i]][9]$beta$psi[,2]) == 2) {		# 2 psi params
	modlist[i,12:15] = NA
} else {

if(length(Models[[i]][9]$beta$psi[,1]) == 3) {		# 3 psi params
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiBa"){
		modlist[i,12] = Models[[i]][9]$beta$psi[3,2] # SE of predator effect
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiA.LitterMass_psiA"){
		modlist[i,13] = Models[[i]][9]$beta$psi[3,2] # SE of litter effect
			} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A4_psiA.I(Acari + Formicidae)_psiA"){
		modlist[i,14] = Models[[i]][9]$beta$psi[3,2] # SE of food effect
			} else {
		modlist[i,15] = Models[[i]][9]$beta$psi[3,2] # SE of seasonal effect
		}}}
} else {
	
if(length(Models[[i]][9]$beta$psi[,1]) == 4) {		# 4 psi params
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.LitterMass_psiA"){
	modlist[i,12] = Models[[i]][9]$beta$psi[3,2] # SE of predator effect
	modlist[i,13] = Models[[i]][9]$beta$psi[4,2] # SE of litter effect
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.I(Acari + Formicidae)_psiA"){
	modlist[i,12] = Models[[i]][9]$beta$psi[3,2] # SE of predator effect
	modlist[i,14] = Models[[i]][9]$beta$psi[4,2] # SE of food effect
		} else {
	modlist[i,12] = Models[[i]][9]$beta$psi[3,2] 
	modlist[i,15] = Models[[i]][9]$beta$psi[4,2] # SE of seasonal effect
    }
   }
  }
 }
}	# End SE if-else loops

}	# End model loop

modlist = as.data.frame(modlist)
head(modlist)

colnames(modlist) = c("ModelName","Pars","AIC",
			"Int","SP","INT2","Litter","Arthropods","Season",
			"seINT","seSP","seINT2","seLitter","seArthropods",
			"seSeason")

modlist[,c("Pars","AIC","Int","SP","INT2","Litter","Arthropods","Season",
	"seINT","seSP","seINT2","seLitter","seArthropods","seSeason")] = 
	lapply(modlist[,c("Pars","AIC","Int","SP","INT2","Litter",
	"Arthropods","Season","seINT","seSP","seINT2","seLitter",
	"seArthropods","seSeason")],
  	function(x) as.numeric(levels(x))[x])
head(modlist)

## 2) Calculate model weights  
##	for each model

modlist = modlist[order(modlist$AIC),]		# Sort by AIC

Sample = 872			# Sample size of sites
AICc = modlist$AIC+((2*modlist$Pars)*(modlist$Pars+1)/(Sample-modlist$Pars-1))
				# AICc; AIC corrected for small sample size
DeltaAICc = AICc - AICc[1]	# DeltaAICc
LL = exp(-0.5*DeltaAICc)	# Log-likelihood
Weight =LL/sum(LL)		# Model weight

ModRes = cbind(modlist[,1:3],AICc,DeltaAICc,LL,Weight,modlist[,4:15])


## 3) Calculate model-averaged coefficients estimates,
##	parameter weights, and unconditional standard error for each parameter!

AvgMod = matrix(NA,3,6, 
	dimnames=list(c("Coefficient","Parameter weight","Unconditional SE"),
	c("Intercept","SP","INT2","Litter","Arthropods","Season")))

subsets = list(c("Int","seINT","Weight"),c("SP","seSP","Weight"),
	c("INT2","seINT2","Weight"),c("Litter","seLitter","Weight"),
	c("Arthropods","seArthropods","Weight"),c("Season","seSeason","Weight"))
n = length(subsets)

for(i in 1:n){			# for each 'subset' from 1-6
	sub = ModRes[,subsets[[i]]]
	AvgMod[1,i] = sum(sub[,1]*sub$Weight,na.rm=TRUE)# Coefficient
	AvgMod[2,i] = sum(subset(sub,sub[,1] != "NA")[,3])# Param weight
	se = (sub[,2]^2+(sub[,1]-AvgMod[1,i])^2)^0.5
	se[is.na(se)] = AvgMod[1,i]
	bind = cbind(sub,se)
	SEsub = subset(bind, bind$se < summary(bind$se)[5])# Remove wild SEs
	AvgMod[3,i] = crossprod(SEsub$se, SEsub$Weight)}	# Unconditional SE
(AvgModOPpsi = AvgMod)

# Copy and paste the detection table into Excel for final manipulation
clip = pipe("pbcopy","w")
write.table(AvgModOPpsi, file=clip, sep="\t", row.names=FALSE)
close(clip)
# TABLE 3


#####
##### DETECTION 
##### ESTIMATES
##### 

## Use for-loops to model-average coefficients for DETECTION

## 1) Extract model name, number of parameters, AIC, and
##    parameter estimates and parameter SEs into a matrix

modlist = matrix(, nrow=N, ncol=17)

for (i in 1:N){		# Model loop

modlist[i,1] = Models[[i]][1]$modname
modlist[i,2] = Models[[i]][7]$npar
modlist[i,3] = Models[[i]][8]$aic
modlist[i,4] = Models[[i]][9]$beta$p[1,1]
modlist[i,5] = Models[[i]][9]$beta$p[2,1]

# Use if-else loops to fill in parameter values for each model
if(length(Models[[i]][9]$beta$p[,1]) == 2) {		# 2 p params
	modlist[i,6:10] = NA
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 3) {		# 3 p params
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_rA[1]"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] # rBA effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].LitterMass_pA"){
		modlist[i,8] = Models[[i]][9]$beta$p[3,1] # Litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].I(Acari + Formicidae)_pA"){
		modlist[i,9] = Models[[i]][9]$beta$p[3,1] # Food effect on detection
	} else {		
		modlist[i,10] = Models[[i]][9]$beta$p[3,1] # Seasonal effect on detection
	}}}
	
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 4) {		# 4 p params
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_rBa[1]"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] # rBA effect on detection
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] # rBa effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].LitterMass_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,8] = Models[[i]][9]$beta$p[4,1] # Litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].I(Acari + Formicidae)_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,9] = Models[[i]][9]$beta$p[4,1] # Food effect on detection
	} else {	
		modlist[i,6] = Models[[i]][9]$beta$p[3,1]  
		modlist[i,10] = Models[[i]][9]$beta$p[4,1] # Seasonal effect on detection
	}}}
	
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 5) {		# 5 p params
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].LitterMass_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] # rBA effect on detection
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] # rBa effect on detection
		modlist[i,8] = Models[[i]][9]$beta$p[5,1] # Litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].I(Acari + Formicidae)_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] 
		modlist[i,9] = Models[[i]][9]$beta$p[5,1] # Food effect on detection
	} else {	
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] 
		modlist[i,10] = Models[[i]][9]$beta$p[5,1] # Seasonal effect on detection
	 }
    }
   }
  }
 }
}	# End parameter if-else loops


modlist[i,11] = Models[[i]][9]$beta$p[1,2]
modlist[i,12] = Models[[i]][9]$beta$p[2,2]

# Use if-else loops to fill in SE values for each model
if(length(Models[[i]][9]$beta$p[,1]) == 2) {		# 2 p params
	modlist[i,13:16] = NA
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 3) {		# 3 p params
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_rA[1]"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] # SE rBA effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].LitterMass_pA"){
		modlist[i,15] = Models[[i]][9]$beta$p[3,2] # SE litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].I(Acari + Formicidae)_pA"){
		modlist[i,16] = Models[[i]][9]$beta$p[3,2] # SE food effect on detection
	} else {	
		modlist[i,17] = Models[[i]][9]$beta$p[3,2] # SE seasonal effect on detection
	}}}
	
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 4) {		# 4 p params
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_rBa[1]"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] # SE of rBA effect on detection
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] # SE of rBa effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].LitterMass_pA"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] # 
		modlist[i,15] = Models[[i]][9]$beta$p[4,2] # SE of litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].I(Acari + Formicidae)_pA"){
		modlist[i,14] = Models[[i]][9]$beta$p[3,2] # 
		modlist[i,16] = Models[[i]][9]$beta$p[4,2] # SE of food effect on detection
	} else {
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,17] = Models[[i]][9]$beta$p[4,2] # SE of  - Arthropods
	}}}
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 5) {		# 5 p params
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].LitterMass_pA"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] 
		modlist[i,15] = Models[[i]][9]$beta$p[5,2] # SE of litter effect
	} else {
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].I(Acari + Formicidae)_pA"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] 
		modlist[i,16] = Models[[i]][9]$beta$p[5,2] # SE of food effect				
	} else {
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] 
		modlist[i,17] = Models[[i]][9]$beta$p[5,2] # SE of season effect
	 }
    }
   }
  }
 }
}	# End SE if-else loops

}	# End model loop

modlist = as.data.frame(modlist)
head(modlist)

colnames(modlist) = c("ModelName","Pars","AIC",
			"Int","SP2","INTo","INTd","Litter","Arthropods","Season",
			"seINT","seSP2","seINTo","seINTd","seLitter","seArthropods","seSeason")

modlist[,c("Pars","AIC","Int","SP2","INTo","INTd","Litter","Arthropods","Season",
	"seINT","seSP2","seINTo","seINTd","seLitter","seArthropods","seSeason")] = 
	lapply(modlist[,c("Pars","AIC","Int","SP2","INTo","INTd","Litter",
	"Arthropods","Season","seINT","seSP2","seINTo","seINTd","seLitter",
	"seArthropods","seSeason")],
  	function(x) as.numeric(levels(x))[x])
head(modlist)

## 2) Calculate model weights  
##	for each model

modlist = modlist[order(modlist$AIC),]		# Sort by AIC

Sample = 872			# Sample size of sites
AICc = modlist$AIC+((2*modlist$Pars)*(modlist$Pars+1)/(Sample-modlist$Pars-1))
				# AICc; AIC corrected for small sample size
DeltaAICc = AICc - AICc[1]	# DeltaAICc
LL = exp(-0.5*DeltaAICc)	# Log-likelihood
Weight = LL/sum(LL)		# Model weight

ModRes = cbind(modlist[,1:3],AICc,DeltaAICc,LL,Weight,modlist[,4:17])
head(ModRes)

## 3) Calculate model-averaged coefficients estimates,
##	parameter weights, and unconditional standard error for each parameter!

AvgMod = matrix(NA,3,7, 
	dimnames=list(c("Coefficient","Parameter weight","Unconditional SE"),
	c("Intercept","SP2","INTo","INTd","Litter","Arthropods","Season")))

subsets = list(c("Int","seINT","Weight"),c("SP2","seSP2","Weight"),
	c("INTo","seINTo","Weight"),c("INTd","seINTd","Weight"),
	c("Litter","seLitter","Weight"),c("Arthropods","seArthropods","Weight"),
	c("Season","seSeason","Weight"))
n = length(subsets)

for(i in 1:n){
	sub = ModRes[,subsets[[i]]]
	AvgMod[1,i] = sum(sub[,1]*sub$Weight,na.rm=TRUE)	# Coefficient
	AvgMod[2,i] = sum(subset(sub,sub[,1] != "NA")[,3]) # Param weight
	se = (sub[,2]^2+(sub[,1]-AvgMod[1,i])^2)^0.5
	se[is.na(se)] = AvgMod[1,i]
	AvgMod[3,i] = crossprod(se,sub$Weight)}				# Unconditional SE
(AvgModOPp = AvgMod)

## Here I generated model-averaged DETECTION parameter estimates,
## parameter weights, and unconditional SEs.

## This information can be used to create:
## 1) A Table ('AvgMod') for the paper with the average model and param weights
## 2) Figures to show how detection varies conditionally, by covariates, or BOTH 

# Copy and paste the detection table into Excel for final manipulation
clip = pipe("pbcopy","w")
write.table(AvgModOPp, file=clip, sep="\t", row.names=FALSE)
close(clip)
# TABLE 4



########################
####### C) NORHUM ######
########################

Models = ModelsNH

####
#### OCCUPANCY
#### ESTIMATES
#### 

## Use for-loops to:
## model-average the coefficients and standard errors

## 1) Extract model name, number of parameters, AIC, and
##    parameter estimates and parameter SEs into a matrix

modlist = matrix(, nrow=N, ncol=15)

for (i in 1:N){

modlist[i,1] = Models[[i]][1]$modname
modlist[i,2] = Models[[i]][7]$npar
modlist[i,3] = Models[[i]][8]$aic
modlist[i,4] = Models[[i]][9]$beta$psi[1,1]
modlist[i,5] = Models[[i]][9]$beta$psi[2,1]

# Use if-else loops to fill in parameter values for each model
if(length(Models[[i]][9]$beta$psi[,1]) == 2) {		# 2 psi params
	modlist[i,6:9] = NA
} else {

if(length(Models[[i]][9]$beta$psi[,1]) == 3) {		# 3 psi params
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiBa"){
		modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # 3A; predator effect on psiA 
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiA.LitterMass_psiA"){
		modlist[i,7] = Models[[i]][9]$beta$psi[3,1] # 3B; litter effect on psiA
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiA.I(Araneae_+_Isopoda)_psiA"){
		modlist[i,8] = Models[[i]][9]$beta$psi[3,1] # 3C; food effect on psiA
		} else {
		modlist[i,9] = Models[[i]][9]$beta$psi[3,1] # 3D; season effect on psiA
		}}}
} else {

if(length(Models[[i]][9]$beta$psi[,1]) == 4) {		# 4 psi params
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.LitterMass_psiA"){
	modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # Predator effect on psiA
	modlist[i,7] = Models[[i]][9]$beta$psi[4,1] # Litter effect on psiA
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.I(Araneae_+_Isopoda)_psiA"){
	modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # Predator effect on psiA
	modlist[i,8] = Models[[i]][9]$beta$psi[4,1] # Food effect on psiA
		} else {
	modlist[i,6] = Models[[i]][9]$beta$psi[3,1] # Predator effect on psiA
	modlist[i,9] = Models[[i]][9]$beta$psi[4,1] # Seasonal effect on psiA
    }
   }
  }
 }
}	# End parameter if-else loops

modlist[i,10] = Models[[i]][9]$beta$psi[1,2]
modlist[i,11] = Models[[i]][9]$beta$psi[2,2]



# Use if-else loops to fill in SE values for each model
if(length(Models[[i]][9]$beta$psi[,2]) == 2) {		# 2 psi params
	modlist[i,12:15] = NA
} else {

if(length(Models[[i]][9]$beta$psi[,1]) == 3) {		# 3 psi params
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiBa"){
		modlist[i,12] = Models[[i]][9]$beta$psi[3,2] # SE of predator effect
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A3_psiA.LitterMass_psiA"){
		modlist[i,13] = Models[[i]][9]$beta$psi[3,2] # SE of litter effect
			} else {
	if(row.names(Models[[i]][9]$beta$psi)[3] == "A4_psiA.I(Araneae_+_Isopoda)_psiA"){
		modlist[i,14] = Models[[i]][9]$beta$psi[3,2] # SE of food effect
			} else {
		modlist[i,15] = Models[[i]][9]$beta$psi[3,2] # SE of seasonal effect
		}}}
} else {
	
if(length(Models[[i]][9]$beta$psi[,1]) == 4) {		# 4 psi params
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.LitterMass_psiA"){
	modlist[i,12] = Models[[i]][9]$beta$psi[3,2] # SE of predator effect
	modlist[i,13] = Models[[i]][9]$beta$psi[4,2] # SE of litter effect
		} else {
	if(row.names(Models[[i]][9]$beta$psi)[4] == "A4_psiA.I(Araneae_+_Isopoda)_psiA"){
	modlist[i,12] = Models[[i]][9]$beta$psi[3,2] # SE of predator effect
	modlist[i,14] = Models[[i]][9]$beta$psi[4,2] # SE of food effect
		} else {
	modlist[i,12] = Models[[i]][9]$beta$psi[3,2] 
	modlist[i,15] = Models[[i]][9]$beta$psi[4,2] # SE of seasonal effect
    }
   }
  }
 }
}	# End SE if-else loops

}	# End model loop

modlist = as.data.frame(modlist)
head(modlist)

colnames(modlist) = c("ModelName","Pars","AIC",
			"Int","SP","INT2","Litter","Arthropods","Season",
			"seINT","seSP","seINT2","seLitter","seArthropods",
			"seSeason")

modlist[,c("Pars","AIC","Int","SP","INT2","Litter","Arthropods","Season",
	"seINT","seSP","seINT2","seLitter","seArthropods","seSeason")] = 
	lapply(modlist[,c("Pars","AIC","Int","SP","INT2","Litter",
	"Arthropods","Season","seINT","seSP","seINT2","seLitter",
	"seArthropods","seSeason")],
  	function(x) as.numeric(levels(x))[x])
head(modlist)

## 2) Calculate model weights  
##	for each model

modlist = modlist[order(modlist$AIC),]		# Sort by AIC

Sample = 872			# Sample size of sites
AICc = modlist$AIC+((2*modlist$Pars)*(modlist$Pars+1)/(Sample-modlist$Pars-1))
				# AICc; AIC corrected for small sample size
DeltaAICc = AICc - AICc[1]	# DeltaAICc
LL = exp(-0.5*DeltaAICc)	# Log-likelihood
Weight =LL/sum(LL)		# Model weight

ModRes = cbind(modlist[,1:3],AICc,DeltaAICc,LL,Weight,modlist[,4:15])


## 3) Calculate model-averaged coefficients estimates,
##	parameter weights, and unconditional standard error for each parameter!

AvgMod = matrix(NA,3,6, 
	dimnames=list(c("Coefficient","Parameter weight","Unconditional SE"),
	c("Intercept","SP","INT2","Litter","Arthropods","Season")))

subsets = list(c("Int","seINT","Weight"),c("SP","seSP","Weight"),
	c("INT2","seINT2","Weight"),c("Litter","seLitter","Weight"),
	c("Arthropods","seArthropods","Weight"),c("Season","seSeason","Weight"))
n = length(subsets)

for(i in 1:n){			# for each 'subset' from 1-6
	sub = ModRes[,subsets[[i]]]
	AvgMod[1,i] = sum(sub[,1]*sub$Weight,na.rm=TRUE)# Coefficient
	AvgMod[2,i] = sum(subset(sub,sub[,1] != "NA")[,3])# Param weight
	se = (sub[,2]^2+(sub[,1]-AvgMod[1,i])^2)^0.5
	se[is.na(se)] = AvgMod[1,i]
	bind = cbind(sub,se)
	SEsub = subset(bind, bind$se < summary(bind$se)[5])# Remove wild SEs
	AvgMod[3,i] = crossprod(SEsub$se, SEsub$Weight)}	# Unconditional SE
(AvgModNHpsi = AvgMod)


# Copy and paste the detection table into Excel for final manipulation
clip = pipe("pbcopy","w")
write.table(AvgModNHpsi, file=clip, sep="\t", row.names=FALSE)
close(clip)
# TABLE 3


##### 
##### DETECTION 
##### ESTIMATES 
##### 

## Use for-loops to model-average coefficients for DETECTION

## 1) Extract model name, number of parameters, AIC, and
##    parameter estimates and parameter SEs into a matrix

modlist = matrix(, nrow=N, ncol=17)

for (i in 1:N){		# Model loop

modlist[i,1] = Models[[i]][1]$modname
modlist[i,2] = Models[[i]][7]$npar
modlist[i,3] = Models[[i]][8]$aic
modlist[i,4] = Models[[i]][9]$beta$p[1,1]
modlist[i,5] = Models[[i]][9]$beta$p[2,1]

# Use if-else loops to fill in parameter values for each model
if(length(Models[[i]][9]$beta$p[,1]) == 2) {		# 2 p params
	modlist[i,6:10] = NA
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 3) {		# 3 p params
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_rA[1]"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] # rBA effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].LitterMass_pA"){
		modlist[i,8] = Models[[i]][9]$beta$p[3,1] # Litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].I(Araneae_+_Isopoda)_pA"){
		modlist[i,9] = Models[[i]][9]$beta$p[3,1] # Food effect on detection
	} else {		
		modlist[i,10] = Models[[i]][9]$beta$p[3,1] # Seasonal effect on detection
	}}}
	
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 4) {		# 4 p params
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_rBa[1]"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] # rBA effect on detection
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] # rBa effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].LitterMass_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,8] = Models[[i]][9]$beta$p[4,1] # Litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].I(Araneae_+_Isopoda)_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,9] = Models[[i]][9]$beta$p[4,1] # Food effect on detection
	} else {	
		modlist[i,6] = Models[[i]][9]$beta$p[3,1]  
		modlist[i,10] = Models[[i]][9]$beta$p[4,1] # Seasonal effect on detection
	}}}
	
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 5) {		# 5 p params
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].LitterMass_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] # rBA effect on detection
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] # rBa effect on detection
		modlist[i,8] = Models[[i]][9]$beta$p[5,1] # Litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].I(Araneae_+_Isopoda)_pA"){
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] 
		modlist[i,9] = Models[[i]][9]$beta$p[5,1] # Food effect on detection
	} else {	
		modlist[i,6] = Models[[i]][9]$beta$p[3,1] 
		modlist[i,7] = Models[[i]][9]$beta$p[4,1] 
		modlist[i,10] = Models[[i]][9]$beta$p[5,1] # Seasonal effect on detection
	 }
    }
   }
  }
 }
}	# End parameter if-else loops


modlist[i,11] = Models[[i]][9]$beta$p[1,2]
modlist[i,12] = Models[[i]][9]$beta$p[2,2]

# Use if-else loops to fill in SE values for each model
if(length(Models[[i]][9]$beta$p[,1]) == 2) {		# 2 p params
	modlist[i,13:16] = NA
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 3) {		# 3 p params
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_rA[1]"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] # SE rBA effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].LitterMass_pA"){
		modlist[i,15] = Models[[i]][9]$beta$p[3,2] # SE litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[3] == "B3_pA[1].I(Araneae_+_Isopoda)_pA"){
		modlist[i,16] = Models[[i]][9]$beta$p[3,2] # SE food effect on detection
	} else {	
		modlist[i,17] = Models[[i]][9]$beta$p[3,2] # SE seasonal effect on detection
	}}}
	
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 4) {		# 4 p params
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_rBa[1]"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] # SE of rBA effect on detection
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] # SE of rBa effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].LitterMass_pA"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] # 
		modlist[i,15] = Models[[i]][9]$beta$p[4,2] # SE of litter effect on detection
	} else {
	if(row.names(Models[[i]][9]$beta$p)[4] == "B4_pA[1].I(Araneae_+_Isopoda)_pA"){
		modlist[i,14] = Models[[i]][9]$beta$p[3,2] # 
		modlist[i,16] = Models[[i]][9]$beta$p[4,2] # SE of food effect on detection
	} else {
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,17] = Models[[i]][9]$beta$p[4,2] # SE of  - Arthropods
	}}}
} else {
if(length(Models[[i]][9]$beta$p[,1]) == 5) {		# 5 p params
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].LitterMass_pA"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] 
		modlist[i,15] = Models[[i]][9]$beta$p[5,2] # SE of litter effect
	} else {
	if(row.names(Models[[i]][9]$beta$p)[5] == "B5_pA[1].I(Araneae_+_Isopoda)_pA"){
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] 
		modlist[i,16] = Models[[i]][9]$beta$p[5,2] # SE of food effect				
	} else {
		modlist[i,13] = Models[[i]][9]$beta$p[3,2] 
		modlist[i,14] = Models[[i]][9]$beta$p[4,2] 
		modlist[i,17] = Models[[i]][9]$beta$p[5,2] # SE of season effect
	 }
    }
   }
  }
 }
}	# End SE if-else loops

}	# End model loop

modlist = as.data.frame(modlist)
head(modlist)

colnames(modlist) = c("ModelName","Pars","AIC",
			"Int","SP2","INTo","INTd","Litter","Arthropods","Season",
			"seINT","seSP2","seINTo","seINTd","seLitter","seArthropods","seSeason")

modlist[,c("Pars","AIC","Int","SP2","INTo","INTd","Litter","Arthropods","Season",
	"seINT","seSP2","seINTo","seINTd","seLitter","seArthropods","seSeason")] = 
	lapply(modlist[,c("Pars","AIC","Int","SP2","INTo","INTd","Litter",
	"Arthropods","Season","seINT","seSP2","seINTo","seINTd","seLitter",
	"seArthropods","seSeason")],
  	function(x) as.numeric(levels(x))[x])
head(modlist)

## 2) Calculate model weights  
##	for each model

modlist = modlist[order(modlist$AIC),]		# Sort by AIC

Sample = 872			# Sample size of sites
AICc = modlist$AIC+((2*modlist$Pars)*(modlist$Pars+1)/(Sample-modlist$Pars-1))
				# AICc; AIC corrected for small sample size
DeltaAICc = AICc - AICc[1]	# DeltaAICc
LL = exp(-0.5*DeltaAICc)	# Log-likelihood
Weight = LL/sum(LL)		# Model weight

ModRes = cbind(modlist[,1:3],AICc,DeltaAICc,LL,Weight,modlist[,4:17])
head(ModRes)

## 3) Calculate model-averaged coefficients estimates,
##	parameter weights, and unconditional standard error for each parameter!

AvgMod = matrix(NA,3,7, 
	dimnames=list(c("Coefficient","Parameter weight","Unconditional SE"),
	c("Intercept","SP2","INTo","INTd","Litter","Arthropods","Season")))

subsets = list(c("Int","seINT","Weight"),c("SP2","seSP2","Weight"),
	c("INTo","seINTo","Weight"),c("INTd","seINTd","Weight"),
	c("Litter","seLitter","Weight"),c("Arthropods","seArthropods","Weight"),
	c("Season","seSeason","Weight"))
n = length(subsets)

for(i in 1:n){
	sub = ModRes[,subsets[[i]]]
	AvgMod[1,i] = sum(sub[,1]*sub$Weight,na.rm=TRUE)	# Coefficient
	AvgMod[2,i] = sum(subset(sub,sub[,1] != "NA")[,3]) # Param weight
	se = (sub[,2]^2+(sub[,1]-AvgMod[1,i])^2)^0.5
	se[is.na(se)] = AvgMod[1,i]
	AvgMod[3,i] = crossprod(se,sub$Weight)}				# Unconditional SE
(AvgModNHp = AvgMod)

## Here I generated model-averaged DETECTION parameter estimates,
## parameter weights, and unconditional SEs.

## This information can be used to create:
## 1) A Table ('AvgMod') for the paper with the average model and param weights
## 2) Figures to show how detection varies conditionally, by covariates, or BOTH 

# Copy and paste the detection table into Excel for final manipulation
clip = pipe("pbcopy","w")
write.table(AvgModNHp, file=clip, sep="\t", row.names=FALSE)
close(clip)
# TABLE 4



##########################################################################
############# IV) Use average models for each species 	      #############
#############     to plot graphs and visualize effects of	  #############
#############	predators on each species  			          #############
##########################################################################

#######################
###### A) CRABRA ###### 
#######################

AvgModCBpsi
AvgModCBp

#### The top models (weight > 0.05) include the following parameters:
CBtopmods[,1]

# Since LitterMass is the only covariate included in the top set,
# plot how each parameters varies by LitterMass

### OCCUPANCY

litter = seq(0,5,0.01)

betaA = AvgModCBpsi[1,1]	# Intercept - psiA
betaBA = AvgModCBpsi[1,2]	# psiBA
betaBa = AvgModCBpsi[1,3]	# psiBa
betaLL = AvgModCBpsi[1,4]	# psi leaf litter

# Model the variation in psiA, psiBA, and psiBa,
# as they relate to variation in leaf litter mass
# using the logit-link: exp(beta)/(1+exp(beta))

# Occupancy of dominant species A
modA = betaA+betaLL*litter
psiA = exp(modA)/(1+exp(modA))	

# Occupancy of species B, given A is present
modBA = betaA+betaBA+betaLL*litter
psiBA = exp(modBA)/(1+exp(modBA))	

# Occupancy of species B, given A is absent
modBa = betaA+betaBA+betaBa+betaLL*litter
psiBa = exp(modBa)/(1+exp(modBa))

## Bind psiA, psiBA, and psiBa together, and plot them all in one graph
occupancy = cbind(litter,psiA,psiBA,psiBa)

# FIGURE 2A
plot(litter, psiA, type="l", lty=1, lwd=4, xlim=c(0.5,4.5), ylim=c(0,1), 
	axes=FALSE,	cex.lab=1.2, cex.axis=1.5,
	xlab="Log(leaf-litter mass [g])", 
	ylab=expression('Proportion of sites occupied '~psi))
axis(1, c(0,0.5,1.5,2.5,3.5,4.5), lwd=4.5, cex.axis=1.3)
axis(2, c(-0.2,0,0.2,0.4,0.6,0.8,1), lwd=4.5, cex.axis=1.3)
lines(litter, psiBA, lty=2, lwd=4)
lines(litter, psiBa, lty=2, lwd=4, col="darkgray")
text(0.7,0.9, "A", cex=2.5)
legend("bottomright", c(expression(psi^A),expression(psi^BA),
	expression(psi^Ba)), inset=0.03, cex=1.1, box.lwd=3,
	lty=c(1,2,2), col=c("black","black","darkgray"),lwd=c(3,3,3))

# Species interaction factor, given variation in leaf-litter
# FIGURE 3A
SIF = (psiA*psiBA)/(psiA*(psiA*psiBA+(1-psiA)*psiBa)) 
plot(litter, SIF, type="l", lty=1, lwd=4, xlim=c(0.5,4.5), ylim=c(0,2), 
	axes=F, cex.lab=1.2, cex.axis=2,
	xlab="Log(leaf-litter mass [g])", ylab="Species interaction factor")
axis(1, c(0,0.5,1.5,2.5,3.5,4.5), lwd=4, cex.axis=1.3)
axis(2, c(-0.1,0,0.33,0.66,1.0,1.33,1.66,2.00), lwd=4, cex.axis=1.3)
abline(h=1, lwd=4, col="darkgray",lty=3)
text(0.7,1.8, "A", cex=2.5)



#### DETECTION
CBtopmods
AvgModCBp

litter = seq(0,5,0.01)

betaPA = AvgModCBp[1,1]		# Intercept - pA
betaPB = AvgModCBp[1,2]		# pB
betaRBA = AvgModCBp[1,3]	# rBA
betaRBa = AvgModCBp[1,4]	# rBa
betaPLL = AvgModCBp[1,5]	# pLL

detA = betaPA+betaPLL*litter	# Detection of species A, 
pA = exp(detA)/(1+exp(detA))	# given species B is absent

detB = betaPA+betaPB+betaPLL*litter	# Detection of species B,
pB = exp(detB)/(1+exp(detB))		# given species A is absent

detBA = betaPA+betaPB+betaRBA+betaPLL*litter	# Detection of species B,
rBA = exp(detBA)/(1+exp(detBA))			# given both species are present,
								# and species A was detected 

detBa = betaPA+betaPB+betaRBA+		# Detection of species B,
	betaRBa+betaPLL*litter			# given both species are present,
rBa = exp(detBa)/(1+exp(detBa))		# and species A was not detected
 
detection = cbind(pA,pB,rBA,rBa)	# Bind all the detection predictions together

# FIGURE 4A
plot(litter, pA, type="l", lty=1, lwd=4, xlim=c(0.5,4.5), ylim=c(0,1), 
	axes=FALSE,	cex.lab=1.2, cex.axis=2,
	xlab="Log(leaf-litter mass [g])", 
	ylab=expression(Detection~probability~italic((p))))
axis(1, c(0,0.5,1.5,2.5,3.5,4.5), lwd=4, cex.axis=1.3)
axis(2, c(-0.2,0,0.2,0.4,0.6,0.8,1), lwd=4, cex.axis=1.3)
lines(litter, pB, lty=2, lwd=4, col="black")
lines(litter, rBA, lty=2, lwd=4, col="darkgray")
#lines(litter, rBa, lty=2, lwd=4, col="darkgrey")
legend("topright", c(expression(italic(p^A)),expression(italic(p^B)),
	expression(italic(r^BA))), inset=0.05, cex=1.1,
	lty=c(1,2,2), col=c("black","black","darkgray"),
	lwd=c(3,3,3), box.lwd=2)
text(0.7,0.9, "A", cex=2.5)

#OR FIGURE 4A WITHOUT LEAF-LITTER COVARIATE, A SIMPLE BAR PLOT
estimates = c(mean(pA),mean(pB),mean(rBA))
barplot(estimates, ylab="Detection probability", axes=FALSE, 
	ylim=c(0,0.6), border=1, names.arg=c(expression(italic(p^A)),
	expression(italic(p^B)),expression(italic(r^BA))))
axis(2, lwd=3, cex.axis=1.3, c(0,0.1,0.2,0.3,0.4,0.5,0.6))
text(0.4,0.55, "A", cex=2.5)
mean(pA)
mean(pB)
mean(rBA)




# Use the 'deltamethod()' to calculate 95% CI for the SIF
library(msm)

(STDER = deltamethod(~ -((exp(x1)/(1+exp(x1)))-1)*(exp(x1+x2+x3)/(1+exp(x1+x2+x3)))/
  (((exp(x1)/(1+exp(x1)))*(exp(x1+x2)/(1+exp(x1+x2)))+
  	(1-(exp(x1)/(1+exp(x1))))*(exp(x1+x2+x3)/(1+exp(x1+x2+x3))))^2),
	topmod$beta$psi[,1],  
	topmod$beta$psi.VC))

#(STDER = deltamethod(~ -((psiA-1)*psiBa)/
#  ((psiA*psiBa)+((1-psiA)*psiBa)^2),
#	topmod$beta$psi[,1],  
#	topmod$beta$psi.VC))

#(CI = 2*STDER)			# 2*SE = ~95% CI

detach(package:msm)


#######################
###### B) OOPPUM ###### 
#######################

AvgModOPpsi
AvgModOPp

#### The top models (weight > 0.05) include the following parameters:
OPtopmods[,1]

# Since LitterMass is the only covariate included in the top set,
# plot how each parameters varies by LitterMass

#### OCCUPANCY

litter = seq(0,5,0.01)

betaA = AvgModOPpsi[1,1]	# Intercept - psiA
betaBA = AvgModOPpsi[1,2]	# psiBA
betaBa = AvgModOPpsi[1,3]	# psiBa
betaLL = AvgModOPpsi[1,4]	# psi leaf litter

# Model the variation in psiA, psiBA, and psiBa,
# as they relate to variation in leaf litter mass
# using the logit-link: exp(beta)/(1+exp(beta))

modA = betaA+betaLL*litter			# psiA
psiA = exp(modA)/(1+exp(modA))

modBA = betaA+betaBA+betaLL*litter		# psiBA
psiBA = exp(modBA)/(1+exp(modBA))	

modBa = betaA+betaBA+betaBa+betaLL*litter	# psiBa
psiBa = exp(modBa)/(1+exp(modBa))

## Bind predictions together, and plot them all in one graph
occupancy = cbind(litter,psiA,psiBA,psiBa)

# FIGURE 2B
plot(litter, psiA, type="l", lty=1, lwd=4, xlim=c(0.5,4.5), ylim=c(0,1), 
	axes=FALSE,	cex.lab=1.2, cex.axis=1.5,
	xlab="Log(leaf-litter mass [g])", 
	ylab=expression('Proportion of sites occupied '~psi))
axis(1, c(0,0.5,1.5,2.5,3.5,4.5), lwd=4.5, cex.axis=1.3)
axis(2, c(-0.2,0,0.2,0.4,0.6,0.8,1), lwd=4.5, cex.axis=1.3)
lines(litter, psiBA, lty=2, lwd=4)
lines(litter, psiBa, lty=2, lwd=4, col="darkgray")
text(0.7,0.9, "B", cex=2.5)
#legend("bottomright", c(expression(psi^A),expression(psi^BA),
#	expression(psi^Ba)), inset=0.03, cex=1.1, box.lwd=3,
#	lty=c(1,2,2), col=c("black","black","darkgray"),lwd=c(3,3,3))
legend("bottomright", c(expression(psi^A),expression(psi^BA),
	expression(psi^Ba)), inset=0.03, cex=1.1, box.lwd=3,
	lty=c(1,2,2), col=c("black","black","darkgray"),lwd=c(3,3,3))

# SIF plot
# FIGURE 3B
SIF = (psiA*psiBA)/(psiA*(psiA*psiBA+(1-psiA)*psiBa)) 
plot(litter, SIF, type="l", lty=1, lwd=4, xlim=c(0.5,4.5), ylim=c(0,2), 
	axes=F, cex.lab=1.2, cex.axis=2,
	xlab="Log(leaf-litter mass [g])", ylab="Species interaction factor")
axis(1, c(0,0.5,1.5,2.5,3.5,4.5), lwd=4, cex.axis=1.3)
axis(2, c(-0.1,0,0.33,0.66,1.0,1.33,1.66,2.00), lwd=4, cex.axis=1.3)
abline(h=1, lwd=4, col="darkgray",lty=3)
text(0.7,1.8, "B", cex=2.5)


#### DETECTION
OPtopmods
AvgModOPp

litter = seq(0,5,0.01)

betaPA = AvgModOPp[1,1]		# Intercept - pA
betaPB = AvgModOPp[1,2]		# pB
betaRBA = AvgModOPp[1,3]	# rBA
betaRBa = AvgModOPp[1,4]	# rBa
betaPLL = AvgModOPp[1,5]	# pLL

detA = betaPA+betaPLL*litter	# Detection of species A, 
pA = exp(detA)/(1+exp(detA))	# given species B is absent

detB = betaPA+betaPB+betaPLL*litter	# Detection of species B,
pB = exp(detB)/(1+exp(detB))		# given species A is absent

detBA = betaPA+betaPB+betaRBA+betaPLL*litter	# Detection of species B,
rBA = exp(detBA)/(1+exp(detBA))			# given both species are present,
								# and species A was detected 

detBa = betaPA+betaPB+betaRBA+		# Detection of species B,
	betaRBa+betaPLL*litter			# given both species are present,
rBa = exp(detBa)/(1+exp(detBa))		# and species A was not detected
 
detection = cbind(pA,pB,rBA,rBa)	# Bind all the detection predictions together

# FIGURE 4B
plot(litter, pA, type="l", lty=1, lwd=4, xlim=c(0.5,4.5), ylim=c(0,1), 
	axes=FALSE,	cex.lab=1.2, cex.axis=2,
	xlab="Log(leaf-litter mass [g])", 
	ylab=expression(Detection~probability~italic((p))))
axis(1, c(0,0.5,1.5,2.5,3.5,4.5), lwd=4, cex.axis=1.3)
axis(2, c(-0.2,0,0.2,0.4,0.6,0.8,1), lwd=4, cex.axis=1.3)
lines(litter, pB, lty=2, lwd=4, col="black")
lines(litter, rBA, lty=2, lwd=4, col="darkgray")
legend("topright", c(expression(italic(p^A)),expression(italic(p^B)),
	expression(italic(r^BA))), inset=0.05, cex=1.2,
	lty=c(1,2,2), col=c("black","black","darkgray"), box.lwd=2,
	lwd=c(3,3,3))
text(0.8,0.9, "B", cex=2.5)

#OR FIGURE 4B WITHOUT LEAF-LITTER COVARIATE, A SIMPLE BAR PLOT
estimates = c(mean(pA),mean(pB),mean(rBA))
barplot(estimates, ylab="Detection probability", axes=FALSE,
	ylim=c(0,0.5),
	names.arg=c(expression(italic(p^A)),expression(italic(p^B)),
	expression(italic(r^BA))))
axis(2, lwd=3, cex.axis=1.3, c(0,0.1,0.2,0.3,0.4,0.5))
text(0.4,0.45, "B", cex=2.5)
mean(pA)
mean(pB)
mean(rBA)


#######################
###### C) NORHUM ###### 
#######################

AvgModNHpsi
AvgModNHp

#### The top models (weight > 0.05) include the following parameters:
AvgModNHpsi[,1]

# Since LitterMass is the only covariate included in the top set,
# plot how each parameters varies by LitterMass

#### OCCUPANCY

litter = seq(0,5,0.01)

betaA = AvgModNHpsi[1,1]	# Intercept - psiA
betaBA = AvgModNHpsi[1,2]	# psiBA
betaBa = AvgModNHpsi[1,3]	# psiBa
betaLL = AvgModNHpsi[1,4]	# psi leaf litter

# Model the variation in psiA, psiBA, and psiBa,
# as they relate to variation in leaf litter mass
# using the logit-link: exp(beta)/(1+exp(beta))

# Occupancy of dominant species A
modA = betaA+betaLL*litter
psiA = exp(modA)/(1+exp(modA))	

# Occupancy of species B, given A is present
modBA = betaA+betaBA+betaLL*litter
psiBA = exp(modBA)/(1+exp(modBA))	

# Occupancy of species B, given A is absent
modBa = betaA+betaBA+betaBa+betaLL*litter
psiBa = exp(modBa)/(1+exp(modBa))

## Bind psiA, psiBA, and psiBa together, and plot them all in one graph
occupancy = cbind(litter,psiA,psiBA,psiBa)

# FIGURE 2C
plot(litter, psiA, type="l", lty=1, lwd=4, xlim=c(0.5,4.5), ylim=c(0,1), 
	axes=FALSE,	cex.lab=1.2, cex.axis=1.5,
	xlab="Log(leaf-litter mass [g])", 
	ylab=expression('Proportion of sites occupied '~psi))
axis(1, c(0,0.5,1.5,2.5,3.5,4.5), lwd=4.5, cex.axis=1.3)
axis(2, c(-0.2,0,0.2,0.4,0.6,0.8,1), lwd=4.5, cex.axis=1.3)
lines(litter, psiBA, lty=2, lwd=4)
lines(litter, psiBa, lty=2, lwd=4, col="darkgray")
text(0.7,0.9, "C", cex=2.5)
legend("bottomright", c(expression(psi^A),expression(psi^BA),
	expression(psi^Ba)), inset=0.03, cex=1.1, box.lwd=3,
	lty=c(1,2,2), col=c("black","black","darkgray"),lwd=c(3,3,3))


# FIGURE 3C
# Species interaction factor, given variation in leaf-litter
SIF = (psiA*psiBA)/(psiA*(psiA*psiBA+(1-psiA)*psiBa)) 
plot(litter, SIF, type="l", lty=1, lwd=4, xlim=c(0.5,4.5), ylim=c(0,6), 
	axes=F, cex.lab=1.2, cex.axis=2,
	xlab="Log(leaf-litter mass [g])", ylab="Species interaction factor")
axis(1, c(0,0.5,1.5,2.5,3.5,4.5), lwd=4, cex.axis=1.3)
axis(2, c(-1,0,1,2,3,4,5,6), lwd=4, cex.axis=1.3)
abline(h=1, lwd=4, col="darkgray",lty=3)
text(0.7,5.5, "C", cex=2.5)
# NOT INCLUDED IN PAPER, BECAUSE psiBa PARAMETER RECEIVED SUCH LITTLE SUPPORT
# MODEL AVERAGED PARAMETER WEIGHT = 0.31


#### DETECTION
NHtopmods
AvgModNHp

litter = seq(0,5,0.01)

betaPA = AvgModNHp[1,1]		# Intercept - pA
betaPB = AvgModNHp[1,2]		# pB
betaRBA = AvgModNHp[1,3]	# rBA
betaRBa = AvgModNHp[1,4]	# rBa
betaPLL = AvgModNHp[1,5]	# pLL

detA = betaPA+betaPLL*litter	# Detection of species A, 
pA = exp(detA)/(1+exp(detA))	# given species B is absent

detB = betaPA+betaPB+betaPLL*litter	# Detection of species B,
pB = exp(detB)/(1+exp(detB))		# given species A is absent

detBA = betaPA+betaPB+betaRBA+betaPLL*litter	# Detection of species B,
rBA = exp(detBA)/(1+exp(detBA))			# given both species are present,
								# and species A was detected 

detBa = betaPA+betaPB+betaRBA+		# Detection of species B,
	betaRBa+betaPLL*litter			# given both species are present,
rBa = exp(detBa)/(1+exp(detBa))		# and species A was not detected
 
detection = cbind(pA,pB,rBA,rBa)	# Bind all the detection predictions together

# FIGURE 4C LINE PLOT
plot(litter, pA, type="l", lty=1, lwd=4, xlim=c(0.5,4.5), ylim=c(0,1), 
	axes=FALSE,	cex.lab=1.2, cex.axis=2,
	xlab="Log(leaf-litter mass [g])", 
	ylab=expression(Detection~probability~italic((p))))
axis(1, c(0,0.5,1.5,2.5,3.5,4.5), lwd=4, cex.axis=1.3)
axis(2, c(-0.2,0,0.2,0.4,0.6,0.8,1), lwd=4, cex.axis=1.3)
lines(litter, pB, lty=2, lwd=4, col="black")
lines(litter, rBA, lty=2, lwd=4, col="darkgrey")
lines(litter, rBa, lty=1, lwd=4, col="darkgrey")
legend("topright", c(expression(italic(p^A)),expression(italic(p^B)),
	expression(italic(r^BA)),expression(italic(r^Ba))), inset=0.05, cex=1.2,
	lty=c(1,2,2,1), col=c("black","black","darkgray","darkgray"),
	lwd=c(3,3,3,3), box.lwd=2)
text(0.7,0.9, "C", cex=2.5)

#OR FIGURE 4C WITHOUT LEAF-LITTER COVARIATE, A SIMPLE BAR PLOT
estimates = c(mean(pA),mean(pB),mean(rBA),mean(rBa))
barplot(estimates, ylab="Detection probability", axes=FALSE, 
	ylim=c(0,0.12), border=1,
	names.arg=c(expression(italic(p^A)),expression(italic(p^B)),
	expression(italic(r^BA)),expression(italic(r^Ba))))
axis(2, lwd=3, cex.axis=1.3, c(0,0.02,0.04,0.06,0.08,0.10,0.12))
text(0.4,0.11, "C", cex=2.5)
mean(pA)
mean(pB)
mean(rBA)
mean(rBa)













#############################################################################
################################## THE END ##################################
#############################################################################