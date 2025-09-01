getwd()
# this command shows the current working directory 
setwd("/home/dc_genomics_r")
dir()
sessionInfo()
date()
Sys.time()
round(3.14)
?round()
# this line creates the object 'first_value' and assigns it the value '1'
first_value <- 1
human_chr_number <- 23
gene_name <- 'pten'
ensemble_url <- 'http://ftp.ensemblgenomes.org/pub/bacteria/release-39/fasta/bacteria_5_collection/escherichia_coli_b_str_rel606/'
human_diploid_chr_number <- 2 *human_chr_number
gene_name <- 'tp53'
rm(gene_name)
chromosome_name <- 'chr02'
od_600_value <- 0.47
chr_position <- '1001701'
spock <- TRUE
pilot <- "Earhart"
snp_genes <- c("OXTR", "ACTN3", "AR", "OPRM1")
snps <- c("rs53576", "rs1815739", "rs6152", "rs1799971")
snp_chromosomes <- c("3", "11", "X", "6")
snp_positions <- c(8762685, 66560624, 67545785, 154039662)
snps[3]
snps[1:3]
snps[c(1, 3, 4)]
snps[c(1:3, 4)]
snp_genes <- c(snp_genes, "CYP1A1", "APOA5")
snp_genes
snp_genes[-6]
snp_genes <- snp_genes[-6]
snp_genes
snp_genes[6] <- "APOA5"
snp_genes
snp_positions[snp_positions > 100000000]
snp_positions > 100000000
snp_positions[c(FALSE, FALSE, FALSE, TRUE)]
which(snp_positions > 100000000)
snp_marker_cutoff <- 100000000
snp_positions[snp_positions > snp_marker_cutoff]
is.na(snp_genes)
c("ACTN3","APOA5", "actn3") %in% snp_genes
snp_genes <- c("OXTR", "ACTN3", "AR", "OPRM1", "CYP1A1", NA, "APOA5")
snp_genes
snp_genes <- snp_genes[-5]
snp_genes <- c(snp_genes, NA, NA)
snp_genes
combined <- c(snp_genes[1], snps[1], snp_chromosomes[1], snp_positions[1])
combined
typeof(combined)
snp_data <- list(genes = snp_genes,
                 refference_snp = snps,
                 chromosome = snp_chromosomes,
                 position = snp_positions)
str(snp_data)
snp_data$position
snp_data$position[1]
variants <- read.csv("/home/mahe/dc_genomics_r/combined_tidy_vcf.csv")
dim(variants)
summary(variants)
subset <- data.frame(variants[, c(1:3, 6)])
mode(variants)
class(variants)
alt_alleles <- subset$ALT
head(alt_alleles)
alt_alleles == "A"
alt_alleles[alt_alleles == "A"]
snps <- c(alt_alleles[alt_alleles == "A"],
          alt_alleles[alt_alleles=="T"],
          alt_alleles[alt_alleles=="G"],
          alt_alleles[alt_alleles=="C"])
