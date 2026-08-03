#Script para limpiar el dataset
library(phyloseq)
library(readr)


## Este codigo es de Alberto
# --- Load OTUs
fileOTU="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/count_table.ESV.4R.csv"
otu.in=read.csv(fileOTU)
otu.pseq=otu_table(as.matrix(otu.in), taxa_are_rows = TRUE)

# --- Load samples metadata
fileSample="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/samples_properties/samples_metadata_AltSubstrSyntax.4R.tsv"
sample_metadata = import_qiime_sample_data(fileSample)
#Que 

# --- Load taxonomy and clean it 

# I try with the manual upload
fileTaxonomy="/home/ajf/Desktop/CNB/phyloseq_tutorial/marine_particles_source_data/sequence_table.ESV.fasta_RDPclassified.txt"
taxa.in=read.table(fileTaxonomy,sep=";") 
colnames(taxa.in)=c("taxa_id","none","Kingdom","sig_Kingdom","Phylum","sig_Phylum","Class","sig_Class",
                    "Order","sig_Order","Family","sig_Family","Genus","sig_Genus")
taxonomy=subset(taxa.in,select=c("taxa_id","Kingdom","Phylum","Class","Order","Family","Genus"))
rownames(taxonomy)=taxonomy$taxa_id
taxonomy=subset(taxonomy,select=-c(taxa_id))
tax.pseq = tax_table(as.matrix(taxonomy))

##final del codigo de Alberto

#Creo el objeto de phyloseq

phylo_obj <- merge_phyloseq(sample_metadata, tax.pseq, otu.pseq)
phylo_obj


#First we remove seawater

phylo_obj_B <- prune_samples(sample_data(phylo_obj)$Media == "Beads", phylo_obj)

phylo_obj_B

#--------------------------------------------------------------------

#Now we take only 1 substatre

phylo_obj_B_alginate <- prune_samples(sample_data(phylo_obj_B)$Substrate == "Alginate",phylo_obj_B)

phylo_obj_B_alginate

#---------------------------------------------------------------------

#Now we take the samples that have at least 5000 reads

min(sample_sums(phylo_obj_B_alginate))

data <- prune_samples(sample_sums(phylo_obj_B_alginate) > 5000, phylo_obj_B_alginate)

data
min(sample_sums(data))

#We take the otus that dont have at least 100 reads

data <- prune_taxa(taxa_sums(data) > 10,data)
data

#Eliminar los otus que no esten en al menos 10 muestras
#---------------------------------------------------------------------------------------------------------
# Paso 1: Obtener la tabla de OTUs (abundancias)
otu_table_data <- otu_table(data)

# Paso 2: Contar cuántas muestras tienen presencia de cada OTU (valores > 0)
otu_presence <- apply(otu_table_data, 1, function(x) sum(x > 0))
?apply
# Paso 3: Filtrar los OTUs presentes en al menos 10 muestras
otus_to_keep <- names(otu_presence[otu_presence >= 10])

# Paso 4: Usar prune_taxa() para quedarte solo con los OTUs filtrados
data2<- prune_taxa(otus_to_keep, data)
#---------------------------------------------------------------------------------------------------------
data2

otu_data_frame <- as.data.frame(otu_table(data2,taxa_are_rows = TRUE))

View(otu_data_frame)



write.table(otu_data_frame, file = "data_filtered_for_sparcc_1096tax.tsv", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t" )


