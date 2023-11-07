library(tidyverse)
library(Biostrings)
library(BiocManager)
library(seqinr)
library(rentrez)
library(randomForest)  #Erik - added randomForest library

#Primary author: Arvind Srinivas
#Contributor information: Erik Dassoff has updated and contributed to this script

###----PART 0: ADJUSTABLE PARAMETERS
#added a few adjustable parameters, but more could be added
ntree <- 150

####----PART 1: DATA EXPLORATION, FILTERING, AND QUALITY CHECKING---- 

#Searching for hits. This function takes in the order_name and returns the search results for carnivora, artiodactyla, and squamata. The terms are already defined in search_terms, where we only want the CytB gene and for it to be 800 to 1200 base pairs in length (CytB gene is on average 1140bp). Also using web_history to gain a large dataset, meaning that it might take time to process. Search_results performs the search based on search_term in the nuccore database.
perform_cytb_search <- function(order_name, retmax = NULL) {
  search_term <- paste(order_name, "[ORGN] AND CytB[Gene] AND 800:1200[SLEN]")
  search_results <- entrez_search(db = "nuccore", term = search_term, use_history = TRUE)
  
  #Data exploration. 
  cat("Search results for", order_name, ":\n")
  cat("Class:", class(search_results), "\n")
  cat("Types of IDs:", search_results$ids, "\n")
  cat("Count of total number of hits:", search_results$count, "\n")
  return(search_results)
}

#Searching through carnivora.
search_results_carnivora <- perform_cytb_search("Carnivora") 

#Searching through artiodactyla.
search_results_artiodactyla <- perform_cytb_search("Artiodactyla")

#Searching through squamata.
search_results_squamata <- perform_cytb_search("Squamata")

#Fetching sequences. Using entrez_fetch function to retrieve carnivora, artiodactyla, and squamata squences based on the web_history. Fetching them in a fasta file format. This may take some time. 
fetch_sequences_and_explore <- function(search_results) {
  sequences_fetch <- entrez_fetch(db = "nuccore", web_history = search_results$web_history, rettype = "fasta")
  
  #Data exploration: checking the class of fetched sequences.
  cat("Class of fetched sequences:", class(sequences_fetch), "\n")
  return(sequences_fetch)
}

#Carnivora fetched sequences. 
carnivora_sequences_fetch <- fetch_sequences_and_explore(search_results_carnivora)

#Artiodactyla fetched sequences. 
artiodactyla_sequences_fetch <- fetch_sequences_and_explore(search_results_artiodactyla)

#Squamata fetched sequences.
squamata_sequences_fetch <- fetch_sequences_and_explore(search_results_squamata)

#Writing all sets of fasta files to a text editor to observe next quality checking steps. Also checking for any unusual sequence lengths, but there shouldn't be too much to worry as SLEN was already set for 800:1200. 
write(carnivora_sequences_fetch, "cytB_carnivora_fetch.fasta", sep = "\n")

write(artiodactyla_sequences_fetch, "cytB_artiodactyla_fetch.fasta", sep = "\n")

write(squamata_sequences_fetch, "cytB_squamata_fetch.fasta", sep = "\n")

#Read all sets of character data into DNAStringSet.
carnivora_stringset <- readDNAStringSet("cytB_carnivora_fetch.fasta")

artiodactyla_stringset <- readDNAStringSet("cytB_artiodactyla_fetch.fasta")

squamata_stringset <- readDNAStringSet("cytB_squamata_fetch.fasta")

#Function that changes all stringsets from previous lines into individual dataframes. cytB_title contains the sequence names and cytB_sequence is the actual sequence of each name. 
stringset_to_dataframe <- function(stringset, df_name) {
  df <- data.frame(cytB_title = names(stringset), cytB_sequence = paste(stringset))
  
  #Data exploration: checking class.
  cat("Class of", df_name, "data frame:", class(df), "\n")
  return(df)
}

#Carnivora dataframe.
carnivora_df <- stringset_to_dataframe(carnivora_stringset, "carnivora")

#Artiodactyla dataframe.
artiodactyla_df <- stringset_to_dataframe(artiodactyla_stringset, "artiodactyla")

#Squamata dataframe.
squamata_df <- stringset_to_dataframe(squamata_stringset, "squamata")

#Adding a column called NucleotideCount to each individual dataframe. It counts the number of characters in each sequence and puts the each count as observations in the new NucleotideCount column.
add_nucleotide_count <- function(df, stringset) {
  df$NucleotideCount <- sapply(stringset, function(seq) nchar(as.character(seq)))
  return(df)
}

#Adding it to carnivora_df.
carnivora_df <- add_nucleotide_count(carnivora_df, carnivora_stringset)

#Adding it to artiodactyla_df.
artiodactyla_df <- add_nucleotide_count(artiodactyla_df, artiodactyla_stringset)

#Adding it to squamata_df.
squamata_df <- add_nucleotide_count(squamata_df, squamata_stringset)

#Adding columns to all dataframes called "Order" and placing the order names with respect to their sequences. 
carnivora_df$Order <- "Carnivora"
artiodactyla_df$Order <- "Artiodactyla"
squamata_df$Order <- "Squamata"

#Create new combined dataframe by combining the carnivora, artiodactyla, and squamata data frames by rows. 
combined_df <- rbind(carnivora_df, artiodactyla_df, squamata_df)

#Filtering steps in the combined_df to create combined_df1. New dataframe combines individual dataframes, as well as removing N's and other symbols from the sequences. The cap was set low (0.0001) to ensure little to no sequences containing Ns. This may take some time to load.
combined_df1 <- combined_df %>%
  mutate(cytB_sequence = str_remove(cytB_sequence, "^[-N]+")) %>%
  mutate(cytB_sequence = str_remove(cytB_sequence, "[-N]+$")) %>%
  mutate(cytB_sequence = str_remove_all(cytB_sequence, "-+")) %>% 
  filter(str_count(cytB_sequence, "N") <= (0.0001 * str_count(cytB_sequence))) 
view(combined_df1) 

class(combined_df1) #Checking class.
dim(combined_df1) #Checking dimensions. There should be 4 columns but many rows. 
table(combined_df1$Order) #Observing counts of order data. Confirming that sequences with N's were removed. 
summary(nchar(combined_df1$cytB_sequence)) #Getting statistics for number of nucleotides in sequence for the CytB gene for each order. 

#Create a faceted histogram for sequence length frequency. Inputting combined_df1 as the dataframe. 
fill_colors <- c("Artiodactyla" = "purple", "Carnivora" = "blue", "Squamata" = "orange")

ggplot(combined_df1, aes(x = NucleotideCount, fill = Order)) +
  geom_histogram(binwidth = 20, color = "black") +
  scale_fill_manual(values = fill_colors) +
  labs(title = "Distribution of Sequence Lengths Between Artiodactyla, Carnivora, and Squamata", x = "Sequence Length (bp)", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5),  strip.text = element_text(size = 12)) + #Centers title.
  facet_wrap(~Order) #Combines all histograms into one. 

####----MAIN ANALYSIS (PART 2): BUILDING THE CLASSIFIER----

#Calculate sequence features. Making a new dataframe called cytB_df. Changing cytB_df to DNAStringset.
cytB_df <- as.data.frame(combined_df1)
cytB_df$cytB_sequence <- DNAStringSet(cytB_df$cytB_sequence)
view(cytB_df)

#Erik - check
class(cytB_df$cytB_sequence)

#Looking at nucleotide frequencies from the cytB_df
cytB_df <- cbind(cytB_df, as.data.frame(letterFrequency(cytB_df$cytB_sequence, letters = c("A", "C", "G", "T"))))

#Proportions of A, C, T, and G into new columns of cytB_df. 
cytB_df$Aproportion <- (cytB_df$A) / (cytB_df$A + cytB_df$T + cytB_df$C + cytB_df$G)

cytB_df$Tproportion <- (cytB_df$T) / (cytB_df$A + cytB_df$T + cytB_df$C + cytB_df$G)

cytB_df$Gproportion <- (cytB_df$G) / (cytB_df$A + cytB_df$T + cytB_df$C + cytB_df$G)

cytB_df$Cproportion <- (cytB_df$C) / (cytB_df$A + cytB_df$T + cytB_df$C + cytB_df$G)

#Use k-mer of length 4 to get tetranucleotide frequency and add these frequencies as a new column in cytB_df. 
cytB_df <- cbind(cytB_df, as.data.frame
                 (oligonucleotideFrequency(x = cytB_df$cytB_sequence, width = 4, as.prob = TRUE)))

cytB_df$cytB_sequence <- as.character(cytB_df$cytB_sequence)

#Take sample from "Order" column of dataframe. 
sample <- min(table(cytB_df$Order))
sample #Lowest number of sequence counts of all the orders. Comes from carnivora. 

validation_size <- floor(0.3 * sample) #Validation size is 30% of sample size.
validation_size

#Reproducibility.
set.seed(7845)

#Create a validation data set by sampling from each category of "Order".
cytB_dfValidation <- cytB_df %>%
  group_by(Order) %>%
  sample_n(validation_size)

#Reproducibility.
set.seed(759385)

#Create a training dataset by filtering out validation samples and sampling the remaining amount of data.
cytB_dfTraining <- cytB_df %>%
  filter(!cytB_title %in% cytB_dfValidation$cytB_title) %>% #Samples not from Validation data.
  group_by(Order) %>%
  sample_n(ceiling(0.7 * sample)) 

#Checking the distribution of orders in the training dataset.
table(cytB_dfTraining$Order)
ncol(cytB_df) #Finding number of columns.

#Making a randomForest classifier with training data. 
order_classifier <- randomForest(
  x = cytB_dfTraining[, 9:268], #Takes proportions and 4-mer possibilities.
  y = as.factor(cytB_dfTraining$Order),
  ntree = ntree,
  importance = TRUE
)

#Viewing results of classifier. This take some time. 
order_classifier 

#Plot the error rate proximity plot with custom line colors.
custom_colors <- c("black", "purple", "blue", "orange")

plot(order_classifier, main = "Error Rate Proximity Plot Based on Order Classification for Training Data", col = custom_colors, lty = c(1,2,2,2), lwd = c(1,1.5,1.5,1.5))
legend("topright", legend = c("Out-of-Bag Samples", "Artiodactyla", "Carnivora", "Squamata"), col = custom_colors, lty = c(1,2,2,2), cex = 0.7)

plot(order_classifier)

#Erik - using the OOB error rate to select the number of trees to be used in the model. First creating separate data frame from order_classifier with new column containing the number of trees associated with the error rate 
oobError <- as.data.frame(order_classifier$err.rate) %>%
  cbind(c(1:ntree))

colnames(oobError)[5] <- "treeNo"

#Erik - finding the minimum OOB error rate and using this to select the number of trees for the model. There can be the same (minimum) OOB error for multiple trees so the second piped element filters for the smallest number of trees to create the simplest/most computationally efficient model. Then, the smallest tree number at which there is the lowest OOB error is selected.
min_treeNo <- oobError %>%
  filter(OOB == min(OOB)) %>%
  filter(treeNo == min(treeNo)) %>%
  select(treeNo)

min_treeNo <- min_treeNo[1,1] %>%
  print()

#Erik - duplicating previous model, but this time with auto-selected No. of trees from OOB error rate. 
order_classifier2 <- randomForest(
  x = cytB_dfTraining[, 9:268], #Takes proportions and 4-mer possibilities.
  y = as.factor(cytB_dfTraining$Order),
  ntree = min_treeNo,
  importance = TRUE
)

#Erik - determining top 10 sequence features were most important in generating predictions (using auto-selected tree No.)
topFeatures <- order_classifier2$importance[1:10, ] %>%
  print()   #Single nucleotide proportions had the most predictive power

#Check on validation data (adjusted for auto-selected tree No.). Erik - Also converting to a vector for downstream analysis.
predictCytBValidation2 <- as.vector(predict(order_classifier2, cytB_dfValidation[,9:268]))

#Erik - Adding identifying numbers to rows, to carry unique sample identifiers through analysis. This will be used for downstream analysis to identify samples 
names(predictCytBValidation2) <- c(1:length(predictCytBValidation2))

#check
names(predictCytBValidation2)

predictCytBValidation2

predictCytBValidation2 #The squamata entries were omitted because it went beyond the max amount of entries to display. 

#Erik - creating a confusion matrix for validation data (using auto-selected tree No.)
validationConfusionMatrix2 <- table(predictCytBValidation2, cytB_dfValidation$Order) %>%
  print()

#Erik - finding the incorrect predictions and filtering data frame to show what the incorrect predictions are for further interrogation to see if there is a reason why this was incorrectly classified. Also including cytB titles, to maintain unique sample identifiers.
incorrectPredictions <- as.data.frame(cbind(predictCytBValidation2, cytB_dfValidation$Order, cytB_dfValidation$cytB_title)) %>%
  filter(predictCytBValidation2 != cytB_dfValidation$Order) 

colnames(incorrectPredictions)[c(2,3)] <- c("dfValidationOrder", "Title")

#Erik - check that column names were added correctly
names(incorrectPredictions)
incorrectPredictions

#Erik - finding the incorrect prediction. Sub-setting for the number of rows and title column returns the value of the title, which is needed for downstream analysis. Afterwards, sub-setting the original data frame for this(or these) titles to find more information on these sequences. 
incorrectPredictions_Title <- incorrectPredictions[1:nrow(incorrectPredictions),"Title"]
incorrectPredictions_Title

incorrectClassInfo <- cytB_dfValidation %>%
  filter(cytB_title == incorrectPredictions_Title) %>%
  view() #the sequence count is within an expected range

#Erik - Here, choosing to select the nucleotide sequence, which I'll copy/paste into BLAST online from NCBI to see if it was possibly mis-labelled or a contaminant, making it difficult to classify: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome. Other quality checks could be performed, like looking for outliers in terms of k-mer proportions.

incorrectClassInfo %>% 
  select(cytB_sequence) %>%
  view()  #Erik - compare BLAST hits to incorrectPredictions object to see if any hits match the expected species. In the current example, there is a high probability that the species is Sus scrofa, so the mis-classification is likely from something other than a mislabelling error, contamination, or an incorrect gene (BLAST identified the CYBRD1 gene and sequence lengths are reasonable). Sequence quality is one possibility, since although the BLAST tool identified it, Sus scrofa was not the only possible species. 
  
#Erik - comparing to validation data without auto-selected tree No.
predictCytBValidation <- predict(order_classifier, cytB_dfValidation[,9:268])
validationConfusionMatrix <- table(predictCytBValidation, cytB_dfValidation$Order) %>%
  print()   #the result is the exact same, but could be more computationally efficient with larger data with fewer trees being trained.


