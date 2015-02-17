# Andorf et al., Analysis of autogating results of public flow cytometry
# experiments in R: A ragweed allergy trial from ImmPort as case study
###############################################################################


##Functions to read in and process the downloaded FLOCK results

##FLOCK.dir:		path to folder in which the taskDownload_xxxx folders are located
##incl.description:	logical: should the Description of the FLOCK populations (e.g. FSC-SSC-FITC CD14-PE CD23-PP CD3+APC CD19lo)
##					be included in the results data.frame
##sum.dupl:			logicle: the %-values will be summed up when the same population ("Description") was detected more than once for one file
read.FLOCK <- function(FLOCK.dir,incl.description=FALSE,sum.dupl=TRUE){
	##Determine folders of the single FLOCK runs:
	FLOCK.run.dirs <- list.dirs(path=FLOCK.dir,recursive=FALSE)

	##Check if at least one folder was found for the provided path:
	if(length(FLOCK.run.dirs)!=0){##At least one folder was found
		##Read in and combine the FLOCK results for each file for each run:
		flock.results <- try(do.call("rbind",lapply(FLOCK.run.dirs,FUN=read.in.FLOCK.runs,sum.dupl=sum.dupl)),silent=TRUE)
		##If an error was returned, no result data.frame is getting returned
		##and print out the error message
		if(inherits(flock.results, "try-error")){
			print("An error occured. Please make sure that all FLOCK result folders in the main folder contain results for the same marker panel.")
			print("Original error message:")
			print(flock.results[[1]])
		}else{
			##Use the factor levels of Description to determine the populationID across all result files:
			flock.results <- cbind(flock.results,populationId=as.numeric(as.factor(flock.results[,"Description"])))
			##Exclude Description column if !incl.description:
			if(!incl.description)flock.results <- flock.results[,which(colnames(flock.results)!="Description")]
			
			return(flock.results)
		}
	}else{
		if(length(list.dirs(path=FLOCK.dir,recursive=TRUE))==0){
			print(paste("Folder",FLOCK.dir,"not found."))
		}else{
			print(paste("No folder detected in",FLOCK.dir,": no results read in"))
		}
	}
}



read.in.FLOCK.runs <- function(curr.dir,sum.dupl){
	##Read in summary.all.txt file in the current FLOCK run result folder:
	##First, check if the summray.all.txt file doesn't exist)
	if(file.exists(paste(curr.dir,"/summary.all.txt",sep=""))){
		##Skip the first 7 rows that contain the header:
		FLOCK.summary.all <- read.table(file=paste(curr.dir,"/summary.all.txt",sep=""),sep=",",header=FALSE,as.is=TRUE,strip.white=TRUE,fill=TRUE,skip=7,col.names=c("FileId","Name","TaskStatus"))
		
		##Read in the props.txt of each fcs file (each a separate folder) of the current FLOCK run and combine them:
		##Apply returns list of data.frames and they get bind by row:
		props.tmp <- apply(FLOCK.summary.all,1,FUN=read.in.props.txt,path=curr.dir,sum.dupl=sum.dupl)
		##If all files have a Task Status different from completed, props.tmp is NULL
		if(!is.null(props.tmp)){
			flock.props <- try(do.call("rbind",props.tmp),silent=TRUE)
			##If an error was returned, set flock.props to NULL (results of this folder will not be added to the data.frame)
			##and print out the error message
			if(inherits(flock.props, "try-error")){
				print(paste(curr.dir,"will be skipped. Check if all subfolders contain FLOCK runs on the same markers"))
				print("Original error message:")
				print(flock.props[[1]])
				flock.props <- NULL
			}
		}else{
			flock.props <- props.tmp
		}
		
		print(paste(curr.dir,"#files:",length(unique(flock.props[,"name"])),", done"))
		return(flock.props)
	}else{
		print(paste(curr.dir,"doesn't contain a summary.all.txt file and wil be skipped."))
		return(NULL)
	}
	
}


read.in.props.txt <- function(FLOCK.summary.fi,path,sum.dupl){
	if(FLOCK.summary.fi["TaskStatus"]=="Task Status: completed"){
		props.tmp <- read.table(file=paste(path,"/",strsplit(FLOCK.summary.fi["FileId"],"File Id: ")[[1]][2],"/props.txt",sep=""),sep="\t",header=TRUE,as.is=TRUE)
		flock.props.tmp <- data.frame(name=paste(strsplit(FLOCK.summary.fi["Name"],"Name: ")[[1]][2],"fcs",sep="."),props.tmp[,colnames(props.tmp)[3:ncol(props.tmp)]])
		##the %values will be summed up when the same FLOCK population ("Description") was detected more than once for this file
		if(sum.dupl){
			if(sum(duplicated(flock.props.tmp[,"Description"]))!=0){
				dupl.description <- unique(flock.props.tmp[duplicated(flock.props.tmp[,"Description"]),"Description"])
				for(dd in 1:length(dupl.description)){
					##Determine the rows with duplicated Description
					dupl.rows <- which(flock.props.tmp[,"Description"]==dupl.description[dd])
					##Replace the percentage of the first detection of the Description with the sum
					flock.props.tmp[dupl.rows[1],"Percentage"] <- sum(flock.props.tmp[dupl.rows,"Percentage"])
				}
				##Exclude the rows of the duplicated detections of the Description
				flock.props.tmp <- flock.props.tmp[!duplicated(flock.props.tmp[,"Description"]),]
			}
		}
		return(flock.props.tmp)
	}else{
		print(paste("Check the task status of file",strsplit(FLOCK.summary.fi["Name"],"Name: ")[[1]][2],":",FLOCK.summary.fi["TaskStatus"]))
		return(NULL)
	}
}







