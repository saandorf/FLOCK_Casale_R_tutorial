# Andorf et al., Analysis of FLOCK autogating results of publicly available
# flow cytometry experiments: A ragweed allergy trial from ImmPort as case
# study and tutorial
###############################################################################

##Load the functions to read in and process the downloaded FLOCK results
source("Functions.R")
##If an error is given containing "cannot open file 'Functions.R': No such file
##or directory", please make sure that
##(a) this MAIN.R file is in the same directory as the Functions.R file
##(b) If (a) is the case, please check the working directory of the current R session by uncommenting 
##    (remove # at beginning of row) and executing the following command:
#getwd()
##    If this path is not to where the two .R files are located, please set the working directory
##    to where the .R files are like in the following row:
#setwd("/User/.../MyRproject")

##Define and create output directory
output.dir <- paste(getwd(),"/Outputs/",sep="")
dir.create(output.dir,showWarnings=FALSE)



#####################################
## Read in FLOCK results
#####################################
##The function read.FLOCK (in Functions.R) is called with the path given to where the FLOCK result folders
##are located.
flock.results <- read.FLOCK(FLOCK.dir=paste(getwd(),"/TutorialData",sep=""))
##If an error occurs, please make sure that the path (FLOCK.dir=) is correct to where your FLOCK result
##folders are located. Also, make sure to unzip the folders if you use the function for your own data. 

##Please note that the provided FLOCK result folders contain only the pops.txt file for each of the raw fcs
##files that were analyzed in FLOCK. For storage space reasons, the other three FLOCK result files
##(flock_results.txt, MFI.txt, population_center.txt) are not provided.
##However, the function works the same way if these files are present in your own data.



#####################################
## Sample Annotation
#####################################
##First, choose which annotation method should be used while the default is to use the TSV files
##(the method that should be used needs to be uncommented (no # at beginning of row) while the
##other one is commented (# at beginning of row)):

##Using the tab delimitated (TSV) files:
annotation.method <- "tab.files"
##Or a ImmPort MySQL database installation when available:
#annotation.method <- "DB"

##Get the annotation data:
if(annotation.method=="tab.files"){
	##Annotate with Tab Separated Value files:
	##These files were downloaded directly from the ImmPort website without further formatting.
	
	##Set path to files and read them in:
	tab.study.files <- paste(getwd(),"/TutorialData/SDY1-DR11_Tab/",sep="")
	file_info <- read.table(file=paste(tab.study.files,"file_info.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)
	expsample_2_file_info <- read.table(file=paste(tab.study.files,"expsample_2_file_info.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)
	biosample_2_expsample <- read.table(file=paste(tab.study.files,"biosample_2_expsample.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)
	biosample <- read.table(file=paste(tab.study.files,"biosample.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)
	arm_2_subject <- read.table(file=paste(tab.study.files,"arm_2_subject.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)
	actual_visit <- read.table(file=paste(tab.study.files,"actual_visit.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)
	planned_visit <- read.table(file=paste(tab.study.files,"planned_visit.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)
	
	##Get vector of unique names of the FCS files from the read in FLOCK results
	unique.fcs.files <- as.character(unique(flock.results[,"name"]))
	
	##Combine the read in files and extract the required metadata for the annotation
	file.annotation.tab <- do.call("rbind",sapply(unique.fcs.files,function(x){
						fi.id <- file_info[which(file_info[,"NAME"]==x),"FILE_INFO_ID"]
						exps <- expsample_2_file_info[expsample_2_file_info[,"FILE_INFO_ID"]==fi.id,"EXPSAMPLE_ACCESSION"]
						bios <- biosample_2_expsample[biosample_2_expsample[,"EXPSAMPLE_ACCESSION"]==exps,"BIOSAMPLE_ACCESSION"]
						subj.av <- biosample[biosample[,"BIOSAMPLE_ACCESSION"]==bios,c("SUBJECT_ACCESSION","ACTUAL_VISIT_ACCESSION")]
						pv <- actual_visit[actual_visit[,"ACTUAL_VISIT_ACCESSION"]==subj.av[,"ACTUAL_VISIT_ACCESSION"],"PLANNED_VISIT_ACCESSION"]
						period <- planned_visit[planned_visit[,"PLANNED_VISIT_ACCESSION"]==pv,"PERIOD_ACCESSION"]
						arm <- arm_2_subject[arm_2_subject[,"SUBJECT_ACCESSION"]==subj.av[,"SUBJECT_ACCESSION"],"ARM_ACCESSION"]
						if(length(period)==0){##No period annotation data available for this file
							return(NULL)
						}else{##Full annotation data available for this file
							return(c(arm_accession=arm,planned_visit_accession=pv,period_accession=period))
						}	
					}))
	##Transform the matrix into a data.frame and 
	##make the file names not be the row names but the entries in the first column called "name"
	file.annotation <- data.frame(cbind(name=rownames(file.annotation.tab),file.annotation.tab),row.names=NULL)
	
	##Get arm names and period titles to be used in the figure legends:
	##Get names for the arm_accessions:
	arm.names <- read.table(file=paste(tab.study.files,"arm_or_cohort.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)[,c("ARM_ACCESSION","NAME")]
	##Get titles for the period_accessions:
	period.title <- read.table(file=paste(tab.study.files,"period.txt",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)[,c("PERIOD_ACCESSION","TITLE")]
	##Get submatrix of period.title with only periods listed for which data is used in this analysis
	period.title <- period.title[period.title[,"PERIOD_ACCESSION"]%in%unique(file.annotation[,"period_accession"]),]
	
	##Get further information for the planned_visit_accessions:
	##Get submatrix of planned_visit with only visits listed for which data is used in this analysis
	pv.info <- planned_visit[planned_visit[,"PLANNED_VISIT_ACCESSION"]%in%unique(file.annotation[,"planned_visit_accession"]),]
	##Extract information about Visit and Week from this to be used in the legend later
	pv.info <- cbind(planned_visit_accession=pv.info[,"PLANNED_VISIT_ACCESSION"],data.frame(visit_number=apply(pv.info,1,function(x){
								split.tmp <- strsplit(as.character(x["VISIT_NAME"]),split=".",fixed=TRUE)[[1]]
								strsplit(split.tmp[grep("Visit_",split.tmp)],"_")[[1]][2]
							})))
	
}else if(annotation.method=="DB"){##The code below is executed if annotation.method is set to "DB"
	##Load the packages DBI and RMySQL or install if not present
	if(!(require(DBI)&require(RMySQL))){
		install.packages('RMySQL')
		library(DBI)
		library(RMySQL)
	}
	##Connect and authenticate to a MySQL database
	##connection parameters need to be adjusted to connect to a personal installation of the database
	con <- dbConnect(MySQL(), user="username", password="pw",dbname="database",host="host")
	
  	##Query the database for the desired metadata for the annotation
	file.annotation <- dbGetQuery(con,paste("SELECT fi.name, a2s.arm_accession, av.planned_visit_accession, pv.period_accession
							FROM file_info fi, expsample_2_file_info e2fi, biosample_2_expsample b2e, biosample b, actual_visit av, planned_visit pv, arm_2_subject a2s
							WHERE fi.name IN ('",paste(unique(flock.results[,"name"]),collapse="','"),"')
							AND fi.file_info_id=e2fi.file_info_id
							AND e2fi.experiment_accession=b2e.experiment_accession
							AND e2fi.expsample_accession=b2e.expsample_accession 
							AND b2e.biosample_accession=b.biosample_accession 
							AND b.actual_visit_accession=av.actual_visit_accession
							AND b.subject_accession=av.subject_accession
							AND av.planned_visit_accession=pv.planned_visit_accession
							AND b.subject_accession=a2s.subject_accession;",sep=""))
	
	##Get arm names and period titles to be used in figure legends:
	##et names for the arm_accessions:
	arm.names <- dbGetQuery(con,paste("SELECT arm_accession, name FROM arm_or_cohort
							WHERE arm_accession IN ('",paste(unique(file.annotation[,"arm_accession"]),collapse="','"),"') ;",sep=""))
	##Using the TSV files for annotation, the column names are in upper cases, so
	##to be consistent, make the column names also upper case for the annotation based on the database
	colnames(arm.names) <- toupper(colnames(arm.names))
	##Get titles for the period_accessions:
	period.title <- dbGetQuery(con,paste("SELECT period_accession, title FROM period
							WHERE period_accession IN ('",paste(unique(file.annotation[,"period_accession"]),collapse="','"),"') ;",sep=""))
	##Using the TSV files for annotation, the column names are in upper cases, so
	##to be consistent, make the column names also upper case for the annotation based on the database
	colnames(period.title) <- toupper(colnames(period.title))
	
	##Get further information for the planned_visit_accessions:
	pv.info <- dbGetQuery(con,paste("SELECT planned_visit_accession, visit_name FROM planned_visit
							WHERE planned_visit_accession IN ('",paste(unique(file.annotation[,"planned_visit_accession"]),collapse="','"),"') ;",sep=""))
	##Extract information about Visit and Week from this to be used in the legend later:
	pv.info <- cbind(planned_visit_accession=pv.info[,"planned_visit_accession"],data.frame(visit_number=apply(pv.info,1,function(x){
								split.tmp <- strsplit(x["visit_name"],split=".",fixed=TRUE)[[1]]
								strsplit(split.tmp[grep("Visit_",split.tmp)],"_")[[1]][2]
							})))
	
	##Close connection:
	dbDisconnect(con)

} else{
	print('Set annotation.method either to "DB" or "tab.files"')
}

##Exclude files for which no complete annotation data is available
flock.results <- flock.results[flock.results[,"name"]%in%file.annotation[,"name"],]



#####################################
## Filtering
#####################################
##1. filter step:
##Exclude populations that were detected by FLOCK in fewer than 25% of the files

##Vector of unique populationIds
u.pop_id <- unique(flock.results[,"populationId"])
##Get vector of populationIds that were detected in >= 25% of the FCS files on which FLOCK was run.
pop.2.keep.NAfilter <- u.pop_id[sapply(u.pop_id,function(x){sum(flock.results[,"populationId"]==x)/length(unique(flock.results[,"name"])) >= 0.25})]


##2. filter step:
##Off the populationIds that were dected in >= 25% of the FCS files,
##keep cell populations that have a coefficient of variation (CV) > 40% in any one arm across all visits

##Make list of FCS file names per arm_accession
files.per.arm <- split(file.annotation$name, file.annotation$arm_accession)
##Determine list of FCS files that passed filter 1 and also filter 2:
pop.2.keep <- pop.2.keep.NAfilter[sapply(pop.2.keep.NAfilter,function(x){
					##Determine rows of the flock.results that belong to the current population_id (x)
					##and extract the name and Percentage column of these
					flock.tmp <- flock.results[flock.results[,"populationId"]==x,c("name","Percentage")]
					
					##"Loop" over the list of FCS file names per arm_accession
					##For each arm_accession calculate the CV of the determined sub data.frame (flock.tmp) of the current populationId
					CV.per.arm <- sapply(files.per.arm,function(y){
								##Percentages of the current arm_accession (y) of the sub data.frame flock.tmp
								perc.val <- flock.tmp[flock.tmp[,"name"]%in%y,"Percentage"]
								##Calculate CV of these values and return
								return(sd(perc.val,na.rm=TRUE)/mean(perc.val,na.rm=TRUE))##When only one value per arm: warning is given
							})
					
					#Return TRUE for the curr. populationId if it has at least one arm_accession with a CV >= 0.40
					return(max(CV.per.arm,na.rm=TRUE)>=.4)
				})]
length(pop.2.keep)##20 populations fulfill both filtering criteria and will be used for the analysis

##Get a data.frame containing only the populations that passed both filters
##(populationIDs that are detected in >= 25% or the FCS files and that have a CV of >= 40% for any one population)
flock.res <- flock.results[which(flock.results[,"populationId"]%in%pop.2.keep),]



#####################################
## Averaging FLOCK population percentages
#####################################

##For each planned_visit-arm_accession-combination, determine the mean of the percentage values for each population.
##This will result in only one value per populationIs for each planned_visit-arm_accession-combination.

##Create a matrix with the columns: arm_accession, planned_visit_accession, period_accession and each of the populationIDs.
##And one row for each planned_visit-arm_accession-combination.
##Fill in arm_accession, planned_visit_accession and period_accession.
arm.pv <- data.frame(arm_accession=sort(unique(file.annotation[,"arm_accession"])),
		planned_visit_accession=rep(sort(unique(file.annotation[,"planned_visit_accession"])),each=length(unique(file.annotation[,c("arm_accession")]))))
flock.arm.pv.means <- cbind(arm.pv,matrix(NA,ncol=length(unique(flock.res[,"populationId"])),dimnames=list(NULL,unique(flock.res[,"populationId"]))))

##Fill in the mean percentage values for all planned_visit-arm_accession-combinations for each population_id
flock.arm.pv.means[,3:ncol(flock.arm.pv.means)] <- t(apply(flock.arm.pv.means,1,function(x){
			##Get FCS file names that belong to the arm and planned_visit for the current row of flock.arm.visit.means:
			files.tmp <- file.annotation[which(file.annotation[,"arm_accession"]==x["arm_accession"] & file.annotation[,"planned_visit_accession"]==x["planned_visit_accession"]),"name"]
			##Mean of percentages of these files per population:
			sapply(names(x)[3:length(x)],function(y){mean(flock.res[which(flock.res[,"name"]%in%files.tmp & flock.res[,"populationId"]==y),"Percentage"],na.rm=TRUE)})
		}))



#####################################
## Principal Component Analysis (PCA)
#####################################

##Perform PCA on this data (the populations means for each visit-arm-combination).
##In this function call, the data gets centered and scaled.
pc <- prcomp(flock.arm.pv.means[,3:ncol(flock.arm.pv.means)], center = TRUE, scale = TRUE)
##For help on the used PCA function, uncomment and run the following line of code:
#?prcomp

##Extract the first two principle component scores
##(formed by projecting the entire data set onto the plane created by the first two components)
scores <- pc$x[,1:2]

##Plot the first two principal component scores
##Use a different color for each ARM
pdf(paste(output.dir,"PCA_arm.pdf",sep=""),width=5,height=5)
par(mar=c(4,4,0.1,0.1))
##Load the RColorBrewer package for further color palettes or install if not present
if(!require("RColorBrewer")){
	install.packages("RColorBrewer")
	library(RColorBrewer)
}
colors <- brewer.pal(length(unique(flock.arm.pv.means[,"arm_accession"])),"BrBG")##colorblind safe according to ColorBrewer2.0 http://www.colorbrewer2.org/
names(colors) <- sort(unique(flock.arm.pv.means[,"arm_accession"]))
ind.colors <- sapply(flock.arm.pv.means[,"arm_accession"], function(x) colors[as.character(x)])
plot(scores[,1], scores[,2], xlab = "PC 1", ylab = "PC 2", main= "",cex=2,lwd=2, pch = 20, col=ind.colors)
legend(x="topright",legend=sapply(names(colors),function(x){paste(arm.names[arm.names[,"ARM_ACCESSION"]==x,"NAME"],sep="")}),col=colors ,pch=19,ncol=1,cex=0.6,title="ARM:")##1.1
dev.off()

##Use a different color for each VISIT and include legends for VISIT colors in the figure.
pdf(paste(output.dir,"PCA_visit.pdf",sep=""),width=5,height=5)
par(mar=c(4,4,0.1,0.1))
colors <- brewer.pal(length(unique(flock.arm.pv.means[,"planned_visit_accession"])),"Set1")
names(colors) <- unique(flock.arm.pv.means[,"planned_visit_accession"])[order(unique(flock.arm.pv.means[,"planned_visit_accession"]))]
ind.colors <- sapply(flock.arm.pv.means[,"planned_visit_accession"], function(x) colors[as.character(x)])
plot(scores[,1], scores[,2], xlab = "PC 1", ylab = "PC 2", main= "",cex=2,lwd=2, pch = 20, col=ind.colors)
legend(x="topright",legend=sapply(names(colors), function(x) paste(pv.info[pv.info[,"planned_visit_accession"]==x,"visit_number"]," (",x,")",sep="")),col=colors ,pch=19,ncol=3,cex=0.6,title="Planned visit number (accession):")##1.1
dev.off()

##Use a different color for PERIOD
pdf(paste(output.dir,"PCA_period.pdf",sep=""),width=5,height=5)
par(mar=c(4,4,0.1,0.1))
colors <- brewer.pal(length(unique(period.title[,"PERIOD_ACCESSION"])),"RdYlBu")##   colorblind safe according to ColorBrewer2.0 http://www.colorbrewer2.org/
#Order decreasing to be ordered from screening (P6) to follow-up (P2)
names(colors) <- period.title[,"PERIOD_ACCESSION"][order(period.title[,"PERIOD_ACCESSION"],decreasing=TRUE)]
pv.2.period <- unique(file.annotation[,c("planned_visit_accession","period_accession")])
ind.colors <- sapply(flock.arm.pv.means[,"planned_visit_accession"], function(x) colors[as.character(pv.2.period[pv.2.period[,"planned_visit_accession"]==x,"period_accession"])])
plot(scores[,1], scores[,2], xlab = "PC 1", ylab = "PC 2", main= "",cex=2,lwd=2, pch = 20, col=ind.colors)
legend(x="topright",legend=sapply(names(colors),function(x){paste(x,": ",period.title[period.title[,"PERIOD_ACCESSION"]==x,"TITLE"],sep="")}),col=colors,pch=19,ncol=1,cex=0.6,title="Period:")##1.1
dev.off()



#####################################
## Statistical analysis of clustering results
#####################################
##Add a column containing the periods to the flock.arm.pv.means table:
flock.arm.pv.p.means <- cbind(flock.arm.pv.means,period_accession=sapply(flock.arm.pv.means[,"planned_visit_accession"], function(x) {pv.2.period[pv.2.period[,"planned_visit_accession"]==x,"period_accession"]}))

##Determine the cell population columns
colnames.populations <- colnames(flock.arm.pv.p.means)[which(!colnames(flock.arm.pv.p.means)%in%c("arm_accession","planned_visit_accession","period_accession"))]

##Perform MANOVA
manova.res <- manova(as.matrix(flock.arm.pv.p.means[,colnames.populations]) ~  period_accession+arm_accession+planned_visit_accession,flock.arm.pv.p.means)
summary(manova.res)

