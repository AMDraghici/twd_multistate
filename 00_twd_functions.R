#############################################
    #Custom R Functions for TWD Project#
#############################################

##########################################
    #Part Zero: General Use Functions#
##########################################

## Concatenate Strings inline 
`%+%` <- function(a,b) paste0(a,b)

## Read in vector of libraries 
## If package is not installed then do so
load_packages <- function(libs, install=TRUE){
  
  if(install==TRUE){
    #Read through vector libraries called by user
    for (pkg in libs) {
      #If package is missing from installed.packages then install it
      if (!pkg %in% rownames(installed.packages())) {
        install.packages(pkg)
      }
    }
  }
  
  ## Attach our libraries
  lapply(libs, require, character.only = TRUE)
}

#############################################
   #Part One: Data Cleaning/Manipulation#
#############################################

#####################
##### Load data #####
#####################

## Read Sighting information
twd_read_sightings <- function(data_dir, filename, year_filter, month_filter){
  
  #Read in data with small processing steps
  sighting_data <- read_csv2(file.path(data_dir,filename)) %>% 
    rowid_to_column("Observation") %>% #add observations by row
    mutate(Date = dmy(Date), 
           Year = year(Date),
           Month = month(Date)) %>% 
    filter(Year %in% year_filter, # Keep select years
           Month %in% month_filter) %>% # Keep observations in select months
    filter(Date != "2008-08-05") # Remove media survey
  
  #Return sightings data
  return(sighting_data)
  
}

## Read Block information 
twd_read_blocks <- function(data_dir, filename){
  
  #Read in data with small processing steps
  block_data <- read_csv2(file.path(data_dir,filename)) %>%
    select(Blocks,
           Block_Start = `Start lat.`,
           Block_End = `End lat.`) %>%
    mutate(Block_Length = Block_Start - Block_End)
  
  #Return block data
  return(block_data)
  
}

## Read Survey information
twd_read_surveys <- function(data_dir, filename, year_filter, month_filter){
  
 #Read in data with small processing steps 
 survey_data <- suppressWarnings(read_csv2(file.path(data_dir,filename))) %>% 
    select(Survey, Date,`Start lat.`,`End lat.`) %>% 
    mutate(Date = dmy(Date),
           Year = year(Date),
           Month = month(Date),
           DOY = yday(Date), # Day of experiment (365 calender)
           Survey_Start=pmin(`Start lat.`,`End lat.`), #order left to right
           Survey_End=pmax(`Start lat.`,`End lat.`)) %>%
    filter(Year %in% year_filter, # Keep surveys for select years
           Month %in% month_filter) %>% # Keep surveys for select months
    group_by(Year) %>%
    mutate(DOE = DOY - min(DOY)+1, # Day of experiment (standardized)
           tstep = lead(DOE) - DOE, # Days between surveys
           Survey = row_number()) %>% # Renumber surveys by year
    ungroup() %>%
    select(-`Start lat.`, -`End lat.`) #Drop extra columns
  
  #Return Survey data
  return(survey_data)
  
}

## Summary of Survey Data
twd_summarize_survey <- function(survey){
  #See which days were surveyed on which year
  survey_summary <- survey %>%
    group_by(Year) %>%
    summarize(n=n(),
              Min = min(DOY),
              Max = max(DOY),
              Total_Days = (Max - Min + 1),
              Missed_Days = (Max - Min + 1) -n)
  #Return survey summary
  return(survey_summary)
}

#######################
##### Format data #####
#######################

## Merge Survey/Blocks/Sightings and Process as Capture-Histories 
twd_merge_data <- function(sightings, blocks, survey, block_filter, days_filter){
  
  #Relevant Columns for survey
  survey <- survey %>%  select(Survey,Date,DOE)
  
  #Join up the datasets and clean 
  twd_data <- crossing(sightings,blocks) %>%
    mutate(BlockTmp = (`Sighting lat.` <= Block_Start)) %>%
    group_by(Observation) %>%
    summarize(Blocks = sum(BlockTmp)) %>%
    left_join(sightings, by = "Observation") %>%
    left_join(survey,by = "Date") %>% ## Add survey number
    filter(Blocks %in% block_filter, ## Limit number of blocks 
           DOE <= days_filter,  ## Reduce number of days per year for testing                          
           Blocks != 0) %>% ## Need to look at these
    arrange(Individual, Year, DOE, Time) %>% 
    group_by(Individual, Year, DOE) %>% ## Create capture histories with blocks
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(Individual, Year, DOE,Survey, Blocks) %>%
    arrange(Individual,Year,DOE)
    
    #Return Cleaned Dataset
    return(twd_data)
      
}

## Compute the overlap for surveys on blocks
twd_compute_overlap <- function(survey, blocks, block_filter, days_filter){
  
  #Filter on only selected blocks
  blocks <- blocks %>% filter(Blocks %in% block_filter)
  
  #Join Survey/Block Data and compute the block coverage/overlap percentage
  survey_blocks <- crossing(survey,blocks) %>%
    mutate(Cover=pmin(Survey_End, Block_Start)-pmax(Survey_Start,Block_End),
           Overlap = pmax(0,Cover/Block_Length),
           Covered = ifelse(is.na(Survey),0,(Overlap > 0))) %>%
    filter(DOE <= days_filter) %>%
    select(Date,Year,Survey,Blocks,Covered,Overlap)
  
  #Return Survey Data with covered blocks
  return(survey_blocks)
  
}

## Drop Inconsistent Entries (individuals who couldnt have been seen)
twd_remove_inconsistent <- function(twd_data,survey_blocks){
  
  ## Check for individuals captured in blocks that weren't covered
  check <- twd_data %>%
    left_join(survey_blocks,by=c("Year","Survey","Blocks")) %>%
    arrange(Year,Survey,Individual) %>%
    filter(Overlap == 0) %>%
    select(Individual,Year,Survey) %>%
    mutate(Inconsistent = TRUE)
  
  ## Remove offending observations
  twd_data <- left_join(twd_data,check,by=c("Individual","Year","Survey")) %>%
    filter(is.na(Inconsistent)) %>%
    select(-Inconsistent) %>%
    mutate(ID = as.integer(as.factor(Individual))) # Add numeric ID variable
  
  #Return Data with inconsistent entries filtered out
  return(twd_data)
  
}

#Create Mapping for subset of blocks selected
twd_create_block_map <- function(block_filter){
  
  #Map subset of blocks to index values
  block_map <- data.frame(Blocks = block_filter, Block_id = 0)
  
  #Loop through and assign values
  for(i in 1:length(block_filter)){
    block_id <- block_filter[i]
    block_map$Block_id[i] <- (which(block_filter == block_id))
  }
  
  #Return map
  return(block_map)
}

###############################
##### Format to JAGS data #####
###############################

twd_create_jags_data <- function(data_dir, filenames, filters, Nmax = 100){
  
  #Extract Filenames
  sighting_name <- filenames[1]
  block_name <-  filenames[2]
  survey_name <- filenames[3] 
  old_female_name <- filenames[4]  
  
  #Extract filters
  year_filter <- filters$years
  month_filter <- filters$months
  block_filter <- filters$blocks
  days_filter <- filters$days
  
  #Set up Error Handling
  if(max(diff(year_filter)) != 1){stop("Years must be in sequence and one apart!")}
  
  #Read in Raw Data and do preliminary processing
  sightings <- suppressMessages(twd_read_sightings(data_dir, sighting_name, year_filter, month_filter))
  blocks <- suppressMessages(twd_read_blocks(data_dir, block_name))
  survey <- suppressMessages(twd_read_surveys(data_dir, survey_name, year_filter, month_filter))
  old_females <- suppressMessages(read_delim(file.path(data_dir,old_female_name), ",")) %>% 
    rename(Individual = `Dolphin Code`) %>% 
    mutate(Individual = substr(Individual,1,6))

  for(i in 1:79){
    old_females[i,1] <- substr(old_females[i,1],1,5)
  }
  
  old_females <- old_females %>% 
    pivot_longer(cols=names(old_females[,2:ncol(old_females)]),
                 names_to = "Year",
                 values_to = "OF") %>% 
    mutate(Year = as.double(Year)) %>% 
    arrange(Individual, Year) 
  
  #Map block numbers to index values (must be seq from 1:|B|)
  block_map <- twd_create_block_map(block_filter)
  
  #Compute survey-block overlap
  survey_blocks <- twd_compute_overlap(survey, blocks, block_filter, days_filter) %>% 
    inner_join(block_map,by = "Blocks") %>% 
    mutate(Blocks = Block_id) %>% 
    select(-Block_id)
  
  #Combine capture histories/apply additional filters/remove inconsistent entries
  twd_data <- twd_merge_data(sightings, blocks, survey, block_filter, days_filter) %>% 
    twd_remove_inconsistent(survey_blocks) %>% 
    inner_join(block_map,by = "Blocks") %>% 
    mutate(Blocks = Block_id) %>% 
    select(-Block_id) 
  
  old_females <- old_females %>%
    filter(Individual %in% unique(twd_data$Individual),
           Year %in% year_filter) %>%
    left_join(distinct(select(twd_data,Individual, ID)))
  
  #Number of days per year
  ndays <- twd_summarize_survey(survey)$n
  #Number of individuals in the data
  nind <- length(unique(twd_data$Individual))
  
  #Create Matrix of observed old female groups
  old_female_matrix <- matrix(nrow = Nmax, ncol = length(year_filter))
  
  for(i in 1:nind){
    for(j in 1:length(year_filter)){
      year <- year_filter[j]
      old_female_matrix[i,j] <- old_females %>% 
        filter(ID == i, Year == year) %>% 
        select(OF) %>% 
        distinct() %>% 
        as.numeric()
    }
  }

  #Create array of observed block information
  block_array <- array(dim=c(nind,length(year_filter),max(ndays)))
  
  for(r in 1:nrow(twd_data)){
    block_array[twd_data$ID[r],
                twd_data$Year[r]-year_filter[1] + 1,
                twd_data$Survey[r]] <- twd_data$Blocks[r]
  }
  
  
  #Create array of observed block information
  overlap <- array(dim=c(length(year_filter),max(ndays),length(block_filter)))
  
  for(r in 1:nrow(survey_blocks)){
    overlap[survey_blocks$Year[r]-year_filter[1] + 1,
            survey_blocks$Survey[r],
            survey_blocks$Blocks[r]] <- survey_blocks$Overlap[r]
  }
  
  #Apply data augmentation
  #Nmax =  Maximum number of individual alive in at least one year
  nextra <- Nmax - nind # Number of additional individuals
  exists <- rep(c(1,NA),c(nind,nextra)) #Unknown for additional individuals
  tmp <- array(NA,dim=c(nextra,length(year_filter),max(ndays)))
  block_array <- abind(block_array,tmp,along = 1)
  
  #Create array of capture indicators
  cap_array <- ifelse(is.na(block_array),0,1)
  
  #Create matrix of known survival information
  alive <- (apply(cap_array,c(1,2),sum) > 0) %>%
    ifelse(1,NA)
  
  for(i in 1:nind){
    first <- min(which(alive[i,]==1))
    last <- max(which(alive[i,]==1))
    alive[i,first:last] <- 1
  }
  
  #Create matrix of known recruitment information
  born <- alive
  born[,length(year_filter)] <- 1
  
  for(i in 1:nind){
    first <- min(which(alive[i,]==1))
    born[i,first:length(year_filter)] <- 1
  }
  
  #Create matrix of known mortality
  dead <- 1 - alive
  dead[,1] <- 0
  
  for(i in 1:nind){
    first <- min(which(alive[i,]==1))
    dead[i,1:first] <- 0
  }
  
  #Create matrix of time steps (days between surveys)
  tstep <- matrix(nrow = length(year_filter), ncol = max(ndays))
  
  for(j in 1:length(year_filter)){
    tstep[j,1:ndays[j]] <- filter(survey,Year == year_filter[j]) %>%
      pull(tstep)
  }

  #If not born state is not known then set OF group to NA
  cap_by_year <- matrix(NA,nrow=Nmax,ncol=length(year_filter))
  
  for(i in 1:length(year_filter)){
    cap_by_year[,i] <- 1*(rowSums(cap_array[,i,]) != 0)
  }
  
  
  #If we did not see you we cannot know your state 
  old_female_matrix[!cap_by_year][old_female_matrix[!cap_by_year] != 1] <- NA 
  old_female_matrix[is.na(born)] <- NA
  
  #If we know that you were zero at t then you must be zero at t-1
  for(i in 1:nrow(old_female_matrix)){
    #If there is a zero replace past NA values with zero 
    
    if(is.na(any(old_female_matrix[i,] == 0))){
      next
    } else if(any(old_female_matrix[i,] == 0)){
      old_female_matrix[i,1:max(which(old_female_matrix[i,] == 0))] <- 0
    #If there are only ones and/or NA values skip the row
    } else {
      next
    }
  }
  
  
  #Add Dummy paired survey matrix (used for subsampled model only)
  survey_pair <- matrix(NA, nrow = length(year_filter), ncol = max(ndays)) 
  
  for(i in 1:nrow(survey_pair)){
    survey_pair[i,1:ndays[i]] <- 0
  }
  
  #Gather up JAGS data
  jags_data <- list(Nyears = length(year_filter),
                    Ndays = ndays,
                    tstep = tstep,
                    Tmax = max(tstep, na.rm = TRUE),
                    Nmax = Nmax,
                    exists = rep(c(1,NA),c(nind,nextra)),
                    block = block_array,
                    capture = cap_array,
                    Nblock = length(block_filter),
                    overlap = overlap,
                    dead = dead,
                    born = born,
                    OF = old_female_matrix, 
                    survey_pair = survey_pair)
  
  #Add prior parameters to JAGS data
  jags_data$alpha_psi <- rep(1,length(block_filter))
  
  #Give a jags_list
  return(jags_data)
}

#Subsample Survey Data
twd_subsample_data <- function(jags_data, tmin = 7, day_start, take2=FALSE){
  
  #Start at Day 1 for each year if user doesn't have a preference
  if(missing(day_start)){day_start <- rep(1, nrow(jags_data$tstep))}
  
  #Empty of Years by Survey Occasions
  survey_keep <- matrix(0,nrow = nrow(jags_data$tstep), ncol = max(jags_data$Ndays))
  #Track which surveys are paired and which are apart
  survey_pair <- matrix(0,nrow = nrow(jags_data$tstep), ncol = max(jags_data$Ndays))
  
  #Loop over years
  for(i in 1:nrow(jags_data$tstep)){
    
    #Count how many days have passed since last survey
    counter <- 0 
    #We keep the first survey we start with
    survey_keep[i,day_start[i]] <- 1
    
    #Loop over Survey Days
    for(j in day_start[i]:(jags_data$Ndays[i]-1)){
      
      #Are we doing windowed surveys? If not then use simple iteration
      if(take2 == FALSE){
        #Add number of steps from j to j+1
        counter <- counter + jags_data$tstep[i,j]
        #If number of days that passed is bigger than minimum period (tmin) then keep survey j+1
        survey_keep[i,j+1] <- 0 + 1*(counter >= tmin)
        #If we keep the survey reset the counter at zero
        counter <- counter*(!counter >= tmin)
        #Windowed Survey: If we just sampled at time j then counter is zero
      } else if(counter == 0) { 
        survey_keep[i,j+1] <- 1 #Keep the paired survey
        survey_pair[i,j+1] <- 1*(jags_data$tstep[i,j] < tmin) #if big enough window consider them seperate
        counter <- counter + jags_data$tstep[i,j+1] #update counter
        next
        #Windowed Survey: Currently counting cases
      } else {
        #Update the keep matrix
        survey_keep[i,j+1] <- 0 + 1*(counter >= tmin) 
        #If number of days that passed is bigger than minimum period (tmin) then keep survey j+1
        counter <- counter + jags_data$tstep[i,j+1]
        #If we keep the survey reset the counter at zero (need to call matrix)
        counter <- counter*(1-survey_keep[i,j+1])
      }
    }
  }
  # Update JAGS data
  
  #Drop tstep and Tmax (not in current model just used for subsampling)
  jags_data <- within(jags_data, rm(tstep))
  jags_data <- within(jags_data, rm(Tmax))
  
  #Update Survey Dependent Variables
  jags_data$Ndays <- rowSums(survey_keep) #New total amount of surveys
  
  #Set up new arrays corresponding to subsampled data
  block_array <- array(dim=c(jags_data$Nmax,jags_data$Nyears,max(jags_data$Ndays)))
  capture_array <- array(dim=c(jags_data$Nmax,jags_data$Nyears,max(jags_data$Ndays)))
  overlap_array <- array(0,dim=c(jags_data$Nyears,max(jags_data$Ndays),jags_data$Nblock))
  sp_matrix <- matrix(NA,nrow = jags_data$Nyears, ncol = max(jags_data$Ndays))
  
  #Populate results from full dataset
  for(i in 1:jags_data$Nyears){
    keep_vector <- which(survey_keep[i,] == 1) #which surveys to keep in year i 
    block_array[,i,1:length(keep_vector)] <- jags_data$block[,i,keep_vector]
    capture_array[,i,1:length(keep_vector)] <- jags_data$capture[,i,keep_vector]
    overlap_array[i,1:length(keep_vector),] <- jags_data$overlap[i,keep_vector,]
    sp_matrix[i,1:length(keep_vector)] <- survey_pair[i,keep_vector] #Paired survey
  }

  #Identify ALL empty capture histories and purge from dataset
  
  #Cases where we never see an animal 
  history_keep <- (apply(capture_array,1:2,sum, na.rm=TRUE) %>% rowSums()) > 0
  
  #Remove all the blank entries
  block_array <- block_array[history_keep,,]

  #Reapply data augmentation
  nind <- sum(history_keep) #Number of spotted individuals
  nextra <- jags_data$Nmax - nind  # Number of additional individuals
  
  #Update Existance Vector
  exists <- rep(c(1,NA),c(nind,nextra)) #Unknown for additional individuals
  
  #Update Block Array
  tmp <- array(NA,dim=c(nextra,jags_data$Nyears,max(jags_data$Ndays)))
  block_array <- abind(block_array,tmp,along = 1)
  
  #Construct New Capture Array (we built it in previous steps for the data augmentation step)
  capture_array <- ifelse(is.na(block_array),0,1)
  
  #Create matrix of known survival information
  alive <- (apply(capture_array,c(1,2),sum) > 0) %>%
    ifelse(1,NA)
  
  #Populate using new data
  for(i in 1:nind){
    first <- min(which(alive[i,]==1))
    last <- max(which(alive[i,]==1))
    alive[i,first:last] <- 1
  }
  
  #Create matrix of known recruitment information
  born <- alive
  born[,jags_data$Nyears] <- 1
  
  for(i in 1:nind){
    first <- min(which(alive[i,]==1))
    born[i,first:jags_data$Nyears] <- 1
  }
  
  #Create matrix of known mortality
  dead <- 1 - alive
  dead[,1] <- 0
  
  for(i in 1:nind){
    first <- min(which(alive[i,]==1))
    dead[i,1:first] <- 0
  }
  
  #Drop blocks that are not surveyed 
  #Careful this doesnt mean we drop blocks that we dont see animals in
  #Means that if a block was not visited (overlap = 0 at all samples)..
  #..with this subset it has to be removed
  block_keep <- apply(overlap_array,3,sum) > 0
  overlap_array <- overlap_array[,,block_keep]
  Nblock <- sum(block_keep)
  alpha_psi <- rep(1, Nblock)
  
  #Update Capture to reflect pairing system
  for(i in 1:jags_data$Nyears){
    #Identify paired surveys 
    paired_surveys <- which(sp_matrix[i,]==1)
    #Skip Cases in which pairs fail to appear
    if(length(paired_surveys) == 0){next}
    #Go through pairs and force zero on paired entry 
    for(j in paired_surveys){
      
      #Update Capture Array
      capture_array[,i,j] <- 
        pmax(0,capture_array[,i,j] - capture_array[,i,j-1])
      
      #Update Block Array (cant know where they are if they were not seen)
      unknown_index <- which(block_array[,i,j]*capture_array[,i,j] == 0)
      block_array[unknown_index,i,j] <- NA
    }
  }

  #Update JAGS data list with subsamples
  jags_data$block <- block_array
  jags_data$capture <- capture_array
  jags_data$overlap <- overlap_array
  jags_data$exists <- exists
  jags_data$dead <- dead
  jags_data$born <- born
  jags_data$Nblock <- Nblock
  jags_data$alpha_psi <- alpha_psi
  jags_data$block_keep <- block_keep
  jags_data$survey_pair <- sp_matrix
  
  #Return data out
  return(jags_data)
}


## Generate initial values
twd_generate_init <- function(jags_data,psi = "flat"){
  
  #Extract Key Values
  block_array <- jags_data$block
  nyears <- jags_data$Nyears
  ndays <- jags_data$Ndays
  blocks <- 1:jags_data$Nblock
  nind <- sum(na.omit(jags_data$exists))
  nextra <- jags_data$Nmax - nind
  dead <- jags_data$dead
  born <- jags_data$born
  
  #Set up Block Data
  block_init <- array(dim = dim(block_array))
  
  #Randomly Sample Blocks
  for(j in 1:length(nyears)){
    block_init[,j,1:ndays[j]] <- ifelse(is.na(block_array[,j,1:ndays[j]]),
                                        sample(blocks),NA)
  }
  
  #Initialize State Vector
  if(psi == "flat"){ #even across all blocks per time step
    
    #Construct Matrix
    psi_init <- rep(1/length(blocks),length(blocks))
    
  } else if(psi == "random"){ #random vector
    
    #Initialize Vector
    x <- runif(length(blocks),0,1)
    psik <- x/sum(x)
    psi_init <- psik
    
  } else if(psi == "time"){ #psi can change over time (randomized)
    #Initialize matrix
    nblock <- length(blocks)
    x <- matrix(runif(nblock*nyears,0,1),nrow = nyears, ncol = nblock,byrow = TRUE)
    psi_init <- t(apply(x,1,function(y)(y/sum(y))))
  }
  
  #Gather initial values for jags
  jags_inits <- list(exists = c(rep(NA,nind),rbinom(nextra,1,.5)),
                     block  = block_init, #random block inits
                     dead = ifelse(is.na(dead),0,NA), #replace dead na with 0 (all alive)
                     born = ifelse(is.na(born),1,NA), #replace born na with 1 (all in)
                    # phi = rep(round(runif(1,0.1,0.9),1),nyears - 1), #random phi
                    # p = rep(round(runif(1,0.1,0.9),1),nyears), #random p 
                     psi= psi_init)  #whichever psi setting was used
 
  #Return list of initial values
  return(jags_inits)
}

################
# Overlap Check#
################

#See how much we are actual visting the blocks with our surveys 
twd_check_overlap <- function(jags_data, block, year_index = NULL){
  
  #Auto populate year index
  if(is.null(year_index)){year_index <- 1:jags_data$Nyears}
  
  p1 <- tibble(Coverage = rowSums(jags_data$overlap[year_index,,block],na.rm=T),
            Years = year_index) %>% ggplot() + 
    geom_col(aes(x=as.factor(Years),y=Coverage)) +
    labs(x = "Years",
         y = "Block Coverage", 
         title ="Yearly Coverage for Block #" %+% block) +
    theme(plot.title = element_text(hjust=0.5))
  
  #Return Plot object
  return(p1)
}


#Plot everything on a grid
twd_check_overlap_yrs <- function(jags_data,block_index=NULL,year_index = NULL){
    
    #Auto populate year index
    if(is.null(year_index)){year_index <- 1:jags_data$Nyears}
  
    #Auto populate block index
    if(is.null(block_index)){block_index <- 1:jags_data$Nblock}
  
    #Store plots
    plot_list <- list()
    
    #Call overlap function for all blocks
    for(i in 1:length(block_index)){
      block <- block_index[i]
      plot_list[[i]] <- twd_check_overlap(jags_data,block,year_index)
    }
  
  #Print the plots
  return(do.call("grid.arrange",c(plot_list)))
  
}

#Coverage across entire study period
twd_check_overlap_agg <- function(jags_data,block_index=NULL,year_index =NULL){
  
  #Auto populate year index
  if(is.null(year_index)){year_index <- 1:jags_data$Nyears}
  
  #Auto populate block index
  if(is.null(block_index)){block_index <- 1:jags_data$Nblock}
  
  #Gather Results
  agg_coverage <- matrix(NA,nrow=length(block_index),ncol=2)
  agg_coverage[,1] <- block_index
  for(i in 1:length(block_index)){
   block <- block_index[i]
   agg_coverage[i,2] <- sum(rowSums(jags_data$overlap[year_index,,block],na.rm=T))
  }
  
  #Create Plot
  p1 <- agg_coverage %>% as_tibble() %>% 
    rename(Blocks = V1, Coverage = V2) %>% 
    ggplot() +
    geom_col(aes(y = as.factor(Blocks), x = Coverage)) +
    labs(y = "Block",
          title ="Total Coverage across Blocks") +
    theme(plot.title = element_text(hjust=0.5)) +
    scale_x_continuous(breaks = seq(0,130,10))
   
   #Show the coverage
   return(p1)
   
}

#############################
# Check Block Capture Rates # 
#############################
twd_check_block <- function(jags_data,year_index = NULL){
  
  #Auto populate year index
  if(is.null(year_index)){year_index <- 1:jags_data$Nyears}
  
  #Gather Results
    capture_counts <- table(jags_data$block[,year_index,]) 
  empty_frame <- tibble(Blocks = 1:jags_data$Nblock,
                         Captures = rep(0,jags_data$Nblock))
  #Create Plot
  p1 <- tibble(Blocks = capture_counts %>% rownames() %>% as.integer(),
               Captures = capture_counts %>% as.numeric()) %>% 
    right_join(empty_frame, by = c("Blocks")) %>% 
    replace_na(list(Captures.x = 0, Captures.y = 0)) %>% 
    mutate(Captures = Captures.x + Captures.y) %>% 
    select(Blocks, Captures) %>% 
    ggplot() +
    geom_col(aes(y= factor(Blocks,levels = c(1:jags_data$Nblock)),
                 x=Captures), col="black") +
    labs(y = "Block",
         x= "Captures",
         title ="Total Captures across Blocks") +
    theme(plot.title = element_text(hjust=0.5)) 
    
  return(p1)
  
}

twd_check_block_yrs <- function(jags_data,year_index = NULL){
  
  #Auto populate year index
  if(is.null(year_index)){year_index <- 1:jags_data$Nyears}
  
  #Store plots
  plot_list <- list()
  
  #Call overlap function for all blocks
  for(i in 1:length(year_index)){
    yr <- year_index[i]
    plot_list[[i]] <- twd_check_block(jags_data,yr) + 
      labs(title = "Captures in Year " %+% yr)
  }
  
  #Print the plots
  return(do.call("grid.arrange",c(plot_list)))
}

twd_check_survey <- function(jags_data,block,year_index = NULL){
  
  #Auto populate year index
  if(is.null(year_index)){year_index <- 1:jags_data$Nyears}
  
  p1 <- tibble(Coverage = 
                 colSums(jags_data$overlap[year_index,,block],na.rm=T),
         Survey = 1:ncol(jags_data$overlap[year_index,,block]),
         Block = block) %>% ggplot() + 
    geom_col(aes(y = Coverage, x = as.factor(Survey))) +
    labs(x = "Survey Number")
  return(p1)
}

twd_check_survey_agg <- function(jags_data, block_index = NULL, year_index = NULL){
 
  #Auto populate year index
  if(is.null(year_index)){year_index <- 1:jags_data$Nyears}
  
  #Auto populate block index
  if(is.null(block_index)){block_index <- 1:jags_data$Nblock}
  
  #Store plots
  plot_list <- list()
  
  #Call overlap function for all blocks
  for(i in 1:length(block_index)){
    block <- block_index[i]
    plot_list[[i]] <- twd_check_survey(jags_data,block,year_index) + 
      labs(title = "Survey Coverage in Block " %+% block) +
      ylim(0,10)
  }
  
  #Print the plots
  return(do.call("grid.arrange",c(plot_list)))
}

#############################################
    #Part Two: Model Fitting in JAGS#
#############################################

#Process Jags in Parallel
twd_run_jags_parallel <- function(jags_data, jags_model, jags_params, par_settings, out_dir, save=TRUE, outname = NULL){
  
  #Make sure you don't run the chain for hours and then lose the data
  if(missing(out_dir)){
    if(save == TRUE){
      stop("Specify an output directory using out_dir or set save == FALSE")
    }
  }
  
  #Create Initial Values for each chain using custom initial function
  jags_init <- list()
  
  #First chain uses flat transition initialization then remaining are random
  for(i in 1:par_settings$n.chains){
    
    if(is.null(par_settings$psi_setting)){
      if(i == 1){psi_setting <- "flat"} else {psi_setting <- "random"} 
    } else {
      psi_setting <- par_settings$psi_setting
    }
    
    jags_init[[i]] <- twd_generate_init(jags_data,psi = psi_setting)
  }

  #Start timer
  timer <- proc.time()
  
  #Assign clusters and pass through dclone
  n.cores <- detectCores() - 1
  workers <- pmin(n.cores, par_settings$n.chains)
  
  cat("Setting up " %+% workers %+% " parallel workers ... \n")
  cl <- makePSOCKcluster(workers)
  tmp <- clusterEvalQ(cl, library(dclone))
  
  cat("Initializing graph nodes and adapting chains with " %+% par_settings$n.adapt %+% " iterations ... \n")
  #Construct graph nodes and initialize in arallel
  parJagsModel(cl = cl, 
               name = 'twd', 
               file = jags_model, 
               data = jags_data,
               n.chains = par_settings$n.chains, 
               n.adapt = par_settings$n.adapt,
               inits = jags_init)
  
  #Burn-in each chain in parallel
  cat("Burning-in chains with "%+% par_settings$n.burn %+% " iterations ... \n")
  parUpdate(cl = cl, 
            object = 'twd', 
            n.iter = par_settings$n.burn)
  
  #Sample from distribution in parallel
  cat("Sampling distribution with " %+% par_settings$n.iter %+% " iterations with thinning of " %+% par_settings$n.thin %+%" ... \n ")
  jags_samples <- parCodaSamples(cl = cl, 
                                 model = 'twd', 
                                 variable.names = jags_params, 
                                 n.iter = par_settings$n.iter, 
                                 thin = par_settings$n.thin)
  
  
  #Release clusters
  stopCluster(cl)
  
  #Stop timer and report time
  time.taken <- proc.time() - timer
  cat("MCMC has finished - time statistics are: ...\n")
  print(time.taken)
  
  ## Save output    
  if(save==TRUE){
    cat("Saving samples and initial values to output directory...\n")
    #Time stamp
    tstamp <- strftime(Sys.time(),format="%Y-%m-%d_%H:%M")
    
    #If no name is given just use the timestamp
    if(is.null(outname)){
      outfile_name <- "/jags_samples_" %+% tstamp
      outinit_name <- "/jags_inits_" %+% tstamp
    } else {
      outfile_name <- "/jags_samples_" %+% outname
      outinit_name <- "/jags_inits_" %+% outname
    }
    
    #Save JAGS Samples
    outfile_samples <- out_dir %+% outfile_name %+% ".rds"
    saveRDS(jags_samples,outfile_samples)
    
    #Save Jags Initial Values
    outfile_inits <- out_dir %+% outinit_name %+% ".rds"
    saveRDS(jags_init,outfile_inits)
  }
  
  #Return MCMC Results
  return(jags_samples)
}


#############################################
      #Part Three: MCMC Diagnostics#
#############################################

#Process GGS file (names can get weird)
twd_construct_ggs <- function(jags_samples,parameter){
  
  #Run ggs
  ggs_object <- ggs(jags_samples,parameter)
  
  #Parameter length
  param_len <- ggs_object %>% select(Parameter) %>% distinct() %>% nrow()
  
  #Make sure that we are filtering on the right thing (eg. not both psi[1:15] and p[1:12])
  if(param_len > 1){
    parameter <- parameter %+% "\\["
    ggs_object <- ggs(jags_samples,parameter)
  }
  #Spit out object
  return(ggs_object)
}

#Extract unconditional betas (lets call it gamma)
twd_uncondition_beta <- function(gg_sample){
  
  #Extract Meaningful Parameters
  beta_name <- as.character(unique(gg_sample$Parameter))
  nchains <- max(unique(gg_sample$Chain))
  niter <- max(unique(gg_sample$Iteration))
  
  #Set up Array for vectorized computation (loops are too slow)
  beta_array <- array(dim=c(niter,length(beta_name),nchains))
  gamma_array <- array(dim=c(niter,length(beta_name),nchains))
  gg_list <- list()
  
  #Go through each chain and apply the recursive formula
  for(i in 1:nchains){
    
    #Grab Beta Simulation output
    beta_array[1:niter,1:length(beta_name),i] <- gg_sample %>% 
      filter(Chain == i) %>% 
      select(-Chain) %>%  
      pivot_wider(id_cols = Iteration, names_from = Parameter, values_from = value) %>% 
      select(-Iteration) %>% as.matrix()
    
    #b1 = g1 
    gamma_array[,1,i] <- beta_array[,1,i]
    
    #Initialize the product of compliments
    beta_comp_prod <- 1
    
    #Apply recursive formula
    for(j in 2:length(beta_name)){
      
      #Update cumulative products
      beta_comp_prod <- beta_comp_prod * (1-beta_array[,j-1,i])
      
      #gj = prod_k=1^j-1(1-bk)*bj
      gamma_array[,j,i] <- beta_comp_prod*beta_array[,j,i]
    }
    
    #Return array to ggs object format
    gam_mat <- gamma_array[,,i]
    colnames(gam_mat) <- beta_name
    gg_list[[i]] <- as_tibble(gam_mat) %>% 
      pivot_longer(cols=beta_name,names_to = "Parameter") %>% 
      mutate(Iteration = sort(rep(1:niter,length(beta_name))),
             Chain = i) %>% 
      select("Iteration", "Chain", "Parameter", "value") %>% 
      arrange(Parameter,Chain,Iteration)
  }
  
  ggs_gamma <- bind_rows(gg_list) %>% 
    mutate(Parameter = factor(Parameter,levels=beta_name,labels=beta_name)) %>% 
    arrange(Parameter,Chain,Iteration)
  attributes(ggs_gamma) <- attributes(gg_sample)
  return(ggs_gamma)
}

#Compute Barplot for effective size of all parameters
twd_effective_size <- function(jags_samples, parameter = NULL){
  
  #Compute Effective Sizes
  par_vals <- effectiveSize(jags_samples)
  par_names <- par_vals %>% names()
  
  #Create tibble
  eff_data <- tibble(Parameter = par_names, value = par_vals)
  
  if(!is.null(parameter)){
    eff_data <- eff_data %>% filter(grepl(parameter, Parameter))
  }
  
  #Plot Results
  p1 <- eff_data %>%
    ggplot() +
    geom_col(aes(y = as.factor(Parameter), x = value), color="black") +
    labs(x = "Effective Sample Sizes",
         y= "Parameters",
         title ="Effective Sample Sizes for MCMC Sampler") +
    theme(plot.title = element_text(hjust=0.5)) 
  
  #Return Plot
  return(p1)
}


#############################################
#Part Four: Inference, Summaries and Figures#
#############################################

#Caterpillar Plot
twd_caterpillar <- function(jags_samples, family=NULL, horiz = F, hpd_range=NULL, x=NULL){

  ## Extract specified variables
  if(!is.null(family) & !is.null(x)){
    print("Values for family and x detected - default to X")
    
    #Extract Parameter names and their indicies
    pars <- x$Parameter
    index <- x$x
    
  } else if (!is.null(family)){
    
    #Filter like ggs
    pars <- grep(family,
                 colnames(jags_samples[[1]]),
                 value = TRUE) 
    #Turn into factor
    pars <- factor(pars,levels=pars)
    #Recover index
    index <- 1:length(pars)
    
  }  else if (!is.null(x)){
    
    #Extract Parameter names and their indicies
    pars <- x$Parameter
    index <- x$x 
    
  } else {
    #Need to filter on something
    stop("Specify family or X!")
  }
    

  ## Compute summary
  if(length(pars) == 1){
    
    summ_list <- summary(jags_samples[,levels(pars)])
    summ <- c(summ_list$statistics,summ_list$quantiles) %>%
      t() %>% 
      as_tibble() %>% 
      add_column(index = index, pars = pars)
    
  } else {
    
    summ <- summary(jags_samples[,levels(pars)]) %>%
      do.call(what = cbind) %>%
      as_tibble() %>%
      add_column(index = index,
                 pars = pars)
  }
  
  ## Construct plot
  p1 <- summ %>%
    ggplot(aes(x = as.factor(index), y = Mean)) +
    geom_point(size=3) +
    geom_segment(aes(x=index,
                     xend=index,
                     y=`25%`,
                     yend=`75%`),lwd=2) +
    geom_segment(aes(x=index,
                     xend=index,
                     y=`2.5%`,
                     yend=`97.5%`)) +
    labs(x = "Parameters",
         y = "HPD",
         title = "Caterpillar Plot for MCMC Results") +
    scale_x_discrete(breaks = as.factor(index), labels = levels(pars)) + 
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5)) 
   
    if(!is.null(hpd_range)){p1 <- p1 + ylim(hpd_range[1],hpd_range[2])}
  
    if(horiz == T){p1 <- p1 + coord_flip()}

  ## Return plot
  return(p1)
}