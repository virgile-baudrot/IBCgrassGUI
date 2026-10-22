# Life span of BE2TH species, mapped to the IBC-grass trait MaxAge (years)
#   annual    -> MaxAge = 1   (dies after its first seed release)
#   biennial  -> MaxAge = 2   (dies after seed release in its second year)
#   perennial -> MaxAge = 100 (value of the original IBC-grass PFT files)
# Sources, in order of priority:
#   1. LEDA Traitbase, plant life span (Kleyer et al. 2008)
#      https://uol.de/en/landeco/research/leda/data-files
#      most frequent record per species, ties resolved to the shortest life span
#   2. WCVP, World Checklist of Vascular Plants (Govaerts et al. 2021), lifeform_description
#      https://sftp.kew.org/pub/data-repositories/WCVP/
#      first life form listed (e.g. "annual or biennial" -> annual)
#   3. manual entry for species found in none of them
# Output: Input-EFSA/BE2TH_lifespan.csv (used by 01_create_BE2TH_FieldEdge.R)
library(data.table)

dBE2TH_full <- read.csv("Input-EFSA/Plant traits full dataset_BE2TH.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
SPECIES <- dBE2TH_full$Species
# names as used in the databases
NAME <- sub(" aggr\\.$", "", SPECIES)
NAME[NAME == "Bidens pilosus"] <- "Bidens pilosa"

MAXAGE <- c(annual=1, biennial=2, perennial=100)

# ---------------------------------------------
# 1. LEDA
# ---------------------------------------------
leda_file <- file.path(tempdir(), "plant_life_span.txt")
if(!file.exists(leda_file)){
  download.file("https://uol.de/f/5/inst/biologie/ag/landeco/download/LEDA/Data_files/plant_life_span.txt", leda_file)
}
# 14 lines of SQL query before the header
dLEDA <- read.csv2(leda_file, skip=14, stringsAsFactors=FALSE, quote="", fileEncoding="latin1")
dLEDA$LifeSpan <- c("annuals"="annual",
                    "strict monocarpic bi-annuals and poly-annuals"="biennial",
                    "perennials"="perennial")[dLEDA$gen..plant.life.span]
dLEDA <- dLEDA[!is.na(dLEDA$LifeSpan), ]
leda_lifespan <- function(name){
  # LEDA writes subspecies as "s."
  x <- dLEDA$LifeSpan[dLEDA$SBS.name == sub(" subsp\\. ", " s. ", name)]
  if(length(x) == 0) return(NA)
  n <- table(factor(x, levels=names(MAXAGE)))
  names(n)[which.max(n)] # which.max takes the first, i.e. shortest, in case of tie
}
LEDA <- sapply(NAME, leda_lifespan)

# ---------------------------------------------
# 2. WCVP
# ---------------------------------------------
wcvp_file <- file.path(tempdir(), "wcvp_names.csv")
if(!file.exists(wcvp_file)){
  download.file("https://sftp.kew.org/pub/data-repositories/WCVP/wcvp.zip", file.path(tempdir(), "wcvp.zip"))
  unzip(file.path(tempdir(), "wcvp.zip"), files="wcvp_names.csv", exdir=tempdir())
}
dWCVP <- fread(wcvp_file, sep="|", quote="",
               select=c("plant_name_id", "taxon_rank", "taxon_status", "taxon_name",
                        "lifeform_description", "accepted_plant_name_id"))
dWCVP <- dWCVP[taxon_rank == "Species"]
# synonyms take the life form of their accepted name, accepted names are matched first
dWCVP$lifeform <- dWCVP$lifeform_description[match(dWCVP$accepted_plant_name_id, dWCVP$plant_name_id)]
dWCVP <- dWCVP[order(dWCVP$taxon_status != "Accepted")]
lifeform <- dWCVP$lifeform[match(NAME, dWCVP$taxon_name)]
lifeform_first <- sub("( or |, ).*", "", lifeform)
WCVP <- ifelse(is.na(lifeform) | lifeform == "", NA,
          ifelse(grepl("biennial", lifeform_first), "biennial",
            ifelse(grepl("annual", lifeform_first), "annual", "perennial")))

# ---------------------------------------------
# 3. combine
# ---------------------------------------------
LifeSpan <- ifelse(!is.na(LEDA), LEDA, WCVP)
Source <- ifelse(!is.na(LEDA), "LEDA", ifelse(!is.na(WCVP), "WCVP", NA))
# manual entries (species simulated but absent or without life form in both databases)
MANUAL <- c("Triticosecale rimpaui"="annual",     # triticale, annual cereal
            "Cerastium comatum"="annual",         # Flora Europaea, Flora of Turkey
            "Achillea crithmifolia"="perennial")  # perennial herb, as almost all Achillea
LifeSpan[SPECIES %in% names(MANUAL)] <- MANUAL[SPECIES[SPECIES %in% names(MANUAL)]]
Source[SPECIES %in% names(MANUAL)] <- "manual"

dLifeSpan <- data.frame(
  Species = SPECIES,
  LifeSpan = LifeSpan,
  MaxAge = MAXAGE[LifeSpan],
  Source = Source
)
table(dLifeSpan$LifeSpan, dLifeSpan$Source, useNA="ifany") # check
write.csv(dLifeSpan, file="Input-EFSA/BE2TH_lifespan.csv", row.names=FALSE)
