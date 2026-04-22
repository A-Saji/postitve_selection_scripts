library(tidyverse)
library(ggrepel)

# Load evolutionary signals
library(dplyr)

setwd("~/angeo_C4")

library(Biostrings)
data(BLOSUM62)
blos <- as.matrix(BLOSUM62)

map_msa_to_maize <- function(msa_seq, maize_seq) {
  
  msa_chars   <- strsplit(msa_seq, "")[[1]]
  maize_chars <- strsplit(maize_seq, "")[[1]]
  
  msa_pos   <- seq_along(msa_chars)
  maize_pos <- rep(NA, length(msa_chars))
  
  i <- 1  # msa pointer
  j <- 1  # maize pointer
  
  while (i <= length(msa_chars) && j <= length(maize_chars)) {
    
    if (msa_chars[i] == "-") {
      i <- i + 1
      next
    }
    
    if (msa_chars[i] == maize_chars[j]) {
      maize_pos[i] <- j
      i <- i + 1
      j <- j + 1
    } else {
      # mismatch → advance maize until match
      j <- j + 1
    }
  }
  
  tibble::tibble(
    msa_position   = msa_pos,
    msa_aa         = msa_chars,
    maize_position = maize_pos
  )
}
  
  

df=read.csv("MSA_positions_table_GRMZM2G171179.csv", header = T)
# Compute consensus across each row
df_consensus <- df%>%
  rowwise() %>%
  mutate(
    Consensus = {
      # Extract all residues in this row except Position column
      residues <- c_across(-Position)
      
      # Remove gaps or missing characters
      residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # If all sequences have gap, return "-"
      if (length(residues) == 0) {
        "-"
      } else {
        # Find most frequent residue
        tbl <- table(residues)
        names(tbl)[which.max(tbl)]
      }
    }
  ) %>%
  ungroup()

# View result
head(df_consensus)
colnames(df_consensus)

# C4 species columns
C4_cols <- c("Sevir.2G082700.1.p","Pahal.4G295000.1","Sevir.2G082500.1.p",
             "GRMZM2G171179","Urofu.4G232000.1","Pavag10G072800.1",
             "Urofu.2G078700.1.p","Urofu.4G231900.1.p","Sevir.4G039400.1",
             "ELECO.r07.6BG0464070.1","ELECO.r07.6AG0512190.1")

# C3 species columns
C3_cols <- c("LOC_Os06g09390.1","Bradi1g46690.3","OEL13225.1","Chala.02G062200.1")     


# Dataset 1: Position, Consensus, and C4 species
df_C4 <- df_consensus %>%
  select(Position, Consensus, all_of(C4_cols))

# Dataset 2: Position, Consensus, and C3 species
df_C3 <- df_consensus %>%
  select(Position, Consensus, all_of(C3_cols))


df_C3_clean <- df_C3 %>%
  filter(if_all(all_of(C3_cols), ~ !is.na(.) & . != "-"))

C3_conserved_positions <- df_C3_clean %>%
  # keep only rows where all C3 columns are identical
  filter(apply(select(., all_of(C3_cols)), 1, function(x) length(unique(x)) == 1)) %>%
  pull(Position)

df_C4_conserved_from_C3 <- df_C4 %>%
  filter(Position %in% C3_conserved_positions)

aa_cols <- setdiff(colnames(df_C4_conserved_from_C3), c("Position", "Consensus"))

n_species <- length(aa_cols)

df_consensus_mean <- df_C4_conserved_from_C3 %>%
  rowwise() %>%
  mutate(
    # collect residues for this row
    residues = list(c_across(all_of(aa_cols))),
    
    # count consensus matches
    Consensus_Matches = sum(unlist(residues) == Consensus, na.rm = TRUE),
    
    # count gaps ("-" or "")
    Gaps = sum(unlist(residues) %in% c("-", ""), na.rm = TRUE),
    
    # compute non-gap sequences
    NonGap_Count = n_species - Gaps,
    
    # compute mean consensus score
    Consensus_Mean = ifelse(
      NonGap_Count > 0,
      Consensus_Matches / NonGap_Count,
      NA_real_
    )
  ) %>%
  ungroup()  

df_consensus_mean <- df_consensus_mean %>%
  mutate(
    LowConsensus = Consensus_Mean < 0.4
  )
library(ggplot2)

ggplot(df_consensus_mean, aes(x = Position, y = Consensus_Mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  
  # Normal points
  geom_point(data = subset(df_consensus_mean, !LowConsensus),
             color = "darkblue", size = 2) +
  
  # Highlighted low-consensus positions
  geom_point(data = subset(df_consensus_mean, LowConsensus),
             color = "red", size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Consensus Mean Across Alignment (positions < 0.4 highlighted)",
    x = "Alignment Position",
    y = "Mean Consensus"
  )

ggplot(df_row_scores, aes(x = Position, y = RowSD)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Standard Deviation of Substitution Scores per Alignment Position",
    x = "Alignment Position",
    y = "Row SD (BLOSUM-based)"
  )
aa_cols <- C4_cols
df_row_scores <- df_consensus_mean %>%
  rowwise() %>%
  mutate(
    RowScore = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely from scoring
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Compute substitution scores only for valid AA
      scores <- sapply(valid_residues, function(aa) {
        if (aa == consensus_aa) {
          0
        } else if (aa %in% colnames(blos)) {
          blos[consensus_aa, aa]
        } else {
          NA
        }
      })
      
      # SUM the scores
      sum(scores, na.rm = TRUE)
    },
    
    RowMean = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Number of valid residues (denominator for mean)
      n_valid <- length(valid_residues)
      
      if (n_valid == 0) {
        NA   # no valid AAs at this position
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        sum(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowMean)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "black", size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution BLOSUM62 Score Across Alignment",
    x = "Alignment Position",
    y = "Row BLOSUM62 Substitution Score"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean < 0.4,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_row_scores$Consensus_Mean=1-df_row_scores$Consensus_Mean

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean > 0.75,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )

df_tmp <- df_row_scores %>%
  mutate(
    Low_Frequency = Consensus_Mean > 0.75,
    Low_Score     = RowMean < -0.5
  )
names(df_tmp)
df_tmp2 <- df_tmp
colnames(df_tmp2)[18]=c("Mean substitution frequency")
colnames(df_tmp2)[21]=c("Mean substitution score")
colnames(df_tmp2)[22]=c("Mean substitution standard deviation")

df_long <- df_tmp2 %>%
  select(
    Position,
    `Mean substitution frequency`,
    `Mean substitution score`,
    `Mean substitution standard deviation`,
    Low_Frequency,
    Low_Score
  ) %>%
  pivot_longer(
    cols = c(
      `Mean substitution frequency`,
      `Mean substitution score`,
      `Mean substitution standard deviation`
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = case_when(
      Metric == "Mean substitution frequency"          ~ Low_Frequency,
      Metric == "Mean substitution score"              ~ Low_Score,
      Metric == "Mean substitution standard deviation" ~ FALSE,
      TRUE ~ FALSE
    )
  )

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

msa_seq <- "MCGGAILSDIIPPP--PPR--RVTAGHLWPESKKPRR-AASGR---RG-APV--EQHEQEEDFEADFEEFEVESGESELESE-DEP-KPFAAPRSALARGGL-NTGAAGVDG--PAANSVKRKRKNQFRGIRRRPWGKWAAEIRDPRKGVRVWLGTFNSPEEAARAYDAEARRIRGKKAKVNFPDEVPTAVSQKRRAAGPA--SLKAPKMDVEEEKPIIKLA----VN-----NMTNSNAYHYPAVVGHNIIPEPFMQTQNMPFAPLVNYAA-----LVNLSSD---QGSNSFGCSDFSLENDSRTPDITSVPAPVATLAAVGESVFVQNTAGHAVASPATGNTGVDLAELEPYM-NFLM-DGGSDDSISTLL--SCDGSQ--DVVSNMDLWSFEDMPMSA-GFYX"
maize_seq <- "MCGGAILSDIIPPPPPRRVTAGHLWPESKKPRRAASGRRGAPVEQHEQEEDFEADFEEFEVESGESELESEDEPKPFAAPRSALARGGLNTGAAGVDGPAANSVKRKRKNQFRGIRRRPWGKWAAEIRDPRKGVRVWLGTFNSPEEAARAYDAEARRIRGKKAKVNFPDEVPTAVSQKRRAAGPASLKAPKMDVEEEKPIIKLAVNNMTNSNAYHYPAVVGHNIIPEPFMQTQNMPFAPLVNYAALVNLSSDQGSNSFGCSDFSLENDSRTPDITSVPAPVATLAAVGESVFVQNTAGHAVASPATGNTGVDLAELEPYMNFLMDGGSDDSISTLLSCDGSQDVVSNMDLWSFEDMPMSAGFY*"
mapping_df_2 <- map_msa_to_maize(msa_seq, maize_seq)

domains <- tribble(
  ~Domain, ~Start, ~End,
  "AP2 DNA-binding domain",	111,	170,
  "consensus disorder prediction",	1,	111,
  #  "Zinc finger SWIM-type profile",	590,	626,
)


domain_ranges_msa_2 <- domains %>%
  rowwise() %>%
  mutate(
    msa_start = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == Start)[1] ],
    msa_end   = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == End)[1] ]
  ) %>%
  ungroup()


df_score <- df_long %>% 
  filter(Metric == "Mean substitution score")

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

df_long_filtered <- df_long %>%
  filter(Metric != "Mean substitution standard deviation")

df_sd_fixed <- df_sd %>%
  mutate(Metric = "Mean substitution score")

maize_labels <- mapping_df_2$maize_position
names(maize_labels) <- mapping_df_2$msa_position
df_long_filtered <- df_long_filtered %>%
  left_join(mapping_df_2,
            by = c("Position" = "msa_position"))


EREB160 <- ggplot(
  df_long,
  aes(x = Position, y = Value, fill = Metric)
) +
  
  ## --- DOMAIN REGIONS (background) ---
  geom_rect(
    data = domain_ranges_msa_2,
    aes(
      xmin = msa_start - 0.5,
      xmax = msa_end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = Domain
    ),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  
  ## --- MAIN BARS (THICKER) ---
  geom_col(
    width = 1,
    alpha = 0.6,
    position = "identity"
  ) +
  ## maize positions
  scale_x_continuous(
    breaks = mapping_df_2$msa_position,
    labels = maize_labels
  ) +
  ## --- HIGHLIGHTED BARS ---
  geom_point(
    data = df_long_filtered %>% filter(Highlight),
    aes(x = Position, y = Value),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = df_long_filtered %>% filter(Highlight),
    aes(label = maize_position),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  ## --- COLORS ---
  scale_fill_manual(
    values = c(
      "Mean substitution frequency" = "steelblue",
      "Mean substitution score"     = "darkred",
      "AP2 DNA-binding domain"      = "red",
      "consensus disorder prediction" = "orange"
    )
  ) +
  
  ## --- AXES ---
  coord_cartesian(ylim = c(-2, 2), expand = FALSE) +
  
  ## --- LABELS & THEME ---
  theme_minimal(base_size = 14) +
  labs(
    title = "EREB160",
    x = "Alignment Position",
    y = "Value",
    fill = "Metric / Domain"
  ) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
EREB160

  

#2 bhlh
df=read.csv("MSA_positions_table_GRMZM2G015666_&_GRMZM2G082586.csv", header = T)
# Compute consensus across each row
df_consensus <- df%>%
  rowwise() %>%
  mutate(
    Consensus = {
      # Extract all residues in this row except Position column
      residues <- c_across(-Position)
      
      # Remove gaps or missing characters
      residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # If all sequences have gap, return "-"
      if (length(residues) == 0) {
        "-"
      } else {
        # Find most frequent residue
        tbl <- table(residues)
        names(tbl)[which.max(tbl)]
      }
    }
  ) %>%
  ungroup()

# View result
head(df_consensus)
colnames(df_consensus)

# C4 species columns
C4_cols <- c("Pavag07G168700.1","Sevir.6G193000.1","Urofu.6G174300.1",
             "Pahal.6G239900.1","ELECO.r07.8BG0663710.1","ELECO.r07.8AG0634730.1",
             "ELECO.r07.6BG0491510.1","ELECO.r07.6AG0538290.1","GRMZM2G082586",
             "GRMZM2G015666_P01","Pavag02G210900.1.p","Urofu.2G244100.1.p",
             "Pahal.2G296800.1.p","Sevir.2G248500.1.p")

# C3 species columns
C3_cols <- c("LOC_Os08g37730.1","Bradi3g38510.1.p","Chala.05G056800.1.p",
             "OEL13113.1","Bradi4g32650.1","LOC_Os09g29360.1","OEL27807.1", "Chala.08G054900.1")     


# Dataset 1: Position, Consensus, and C4 species
df_C4 <- df_consensus %>%
  select(Position, Consensus, all_of(C4_cols))

# Dataset 2: Position, Consensus, and C3 species
df_C3 <- df_consensus %>%
  select(Position, Consensus, all_of(C3_cols))


df_C3_clean <- df_C3 %>%
  filter(if_all(all_of(C3_cols), ~ !is.na(.) & . != "-"))

C3_conserved_positions <- df_C3_clean %>%
  # keep only rows where all C3 columns are identical
  filter(apply(select(., all_of(C3_cols)), 1, function(x) length(unique(x)) == 1)) %>%
  pull(Position)

df_C4_conserved_from_C3 <- df_C4 %>%
  filter(Position %in% C3_conserved_positions)

aa_cols <- setdiff(colnames(df_C4_conserved_from_C3), c("Position", "Consensus"))

n_species <- length(aa_cols)

df_consensus_mean <- df_C4_conserved_from_C3 %>%
  rowwise() %>%
  mutate(
    # collect residues for this row
    residues = list(c_across(all_of(aa_cols))),
    
    # count consensus matches
    Consensus_Matches = sum(unlist(residues) == Consensus, na.rm = TRUE),
    
    # count gaps ("-" or "")
    Gaps = sum(unlist(residues) %in% c("-", ""), na.rm = TRUE),
    
    # compute non-gap sequences
    NonGap_Count = n_species - Gaps,
    
    # compute mean consensus score
    Consensus_Mean = ifelse(
      NonGap_Count > 0,
      Consensus_Matches / NonGap_Count,
      NA_real_
    )
  ) %>%
  ungroup()  

df_consensus_mean <- df_consensus_mean %>%
  mutate(
    LowConsensus = Consensus_Mean < 0.4
  )
library(ggplot2)

ggplot(df_consensus_mean, aes(x = Position, y = Consensus_Mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  
  # Normal points
  geom_point(data = subset(df_consensus_mean, !LowConsensus),
             color = "darkblue", size = 2) +
  
  # Highlighted low-consensus positions
  geom_point(data = subset(df_consensus_mean, LowConsensus),
             color = "red", size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Consensus Mean Across Alignment (positions < 0.4 highlighted)",
    x = "Alignment Position",
    y = "Mean Consensus"
  )

ggplot(df_row_scores, aes(x = Position, y = RowSD)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Standard Deviation of Substitution Scores per Alignment Position",
    x = "Alignment Position",
    y = "Row SD (BLOSUM-based)"
  )
aa_cols <- C4_cols
df_row_scores <- df_consensus_mean %>%
  rowwise() %>%
  mutate(
    RowScore = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely from scoring
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Compute substitution scores only for valid AA
      scores <- sapply(valid_residues, function(aa) {
        if (aa == consensus_aa) {
          0
        } else if (aa %in% colnames(blos)) {
          blos[consensus_aa, aa]
        } else {
          NA
        }
      })
      
      # SUM the scores
      sum(scores, na.rm = TRUE)
    },
    
    RowMean = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Number of valid residues (denominator for mean)
      n_valid <- length(valid_residues)
      
      if (n_valid == 0) {
        NA   # no valid AAs at this position
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        sum(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowMean)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "black", size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution BLOSUM62 Score Across Alignment",
    x = "Alignment Position",
    y = "Row BLOSUM62 Substitution Score"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean < 0.4,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_row_scores$Consensus_Mean=1-df_row_scores$Consensus_Mean

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean > 0.75,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )

df_tmp <- df_row_scores %>%
  mutate(
    Low_Frequency = Consensus_Mean > 0.75,
    Low_Score     = RowMean < -0.5
  )
names(df_tmp)
df_tmp2 <- df_tmp
colnames(df_tmp2)[21]=c("Mean substitution frequency")
colnames(df_tmp2)[24]=c("Mean substitution score")
colnames(df_tmp2)[25]=c("Mean substitution standard deviation")

df_long <- df_tmp2 %>%
  select(
    Position,
    `Mean substitution frequency`,
    `Mean substitution score`,
    `Mean substitution standard deviation`,
    Low_Frequency,
    Low_Score
  ) %>%
  pivot_longer(
    cols = c(
      `Mean substitution frequency`,
      `Mean substitution score`,
      `Mean substitution standard deviation`
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = case_when(
      Metric == "Mean substitution frequency"          ~ Low_Frequency,
      Metric == "Mean substitution score"              ~ Low_Score,
      Metric == "Mean substitution standard deviation" ~ FALSE,
      TRUE ~ FALSE
    )
  )

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

msa_seq <- "MALEAVVLSQSQF-QQ--GS-RFG---CG-AM-AA-----G----GAWS-DL-L------FS---G---TE-GLFE-MG--G---A----AGRGWNAA---AS--SPPQ-LL--------QELGD--N-----GAAG-------SLPVASAGS----ASG----------------------------------------GAG---QDAP-------AMAAA--AA-SGRR-KRRRMRPVKNEEEVESQRMIHIAVERNRRKQMNEHLAALRSLMPPAHTQR---------------------------------GDQASIVGGAINFVKELEQLLQSLEARRR-SPQ--C--A-AY---AV-----DP-DDAGPFADFLTFPQYSMCAVIA-APE-NTGH-------------------H-R-EG--G-AVAE--QEAS-GSKPSAVADVEATMVESHANLRVLSRRRPRQLLRLVLGLQGHRLTVLHLNMSS--GAHMVLYSFSLKVEDDCQLTSVGEIAAAAHHIVEKINEEQ-N--------------------K------------A--A----A*"
maize_seq <- "MALEAVVLSQSQFQQGSRFGCGAMAAGGAWSDLLFSGTEGLFEMGGAAGRGWNAAASSPPQLLQELGDNGAAGSLPVASAGSASGGAGQDAPAMAAAAASGRRKRRRMRPVKNEEEVESQRMIHIAVERNRRKQMNEHLAALRSLMPPAHTQRGDQASIVGGAINFVKELEQLLQSLEARRRSPQCAAYAVDPDDAGPFADFLTFPQYSMCAVIAAPENTGHHREGGAVAEQEASGSKPSAVADVEATMVESHANLRVLSRRRPRQLLRLVLGLQGHRLTVLHLNMSSGAHMVLYSFSLKVEDDCQLTSVGEIAAAAHHIVEKINEEQNKAAA*"
mapping_df_2 <- map_msa_to_maize(msa_seq, maize_seq)

domains <- tribble(
  ~Domain, ~Start, ~End,
  "Helix-loop-helix DNA-binding domain",	120,	171,
  "consensus disorder prediction",	80,	116,
  #  "Zinc finger SWIM-type profile",	590,	626,
)


domain_ranges_msa_2 <- domains %>%
  rowwise() %>%
  mutate(
    msa_start = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == Start)[1] ],
    msa_end   = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == End)[1] ]
  ) %>%
  ungroup()


df_score <- df_long %>% 
  filter(Metric == "Mean substitution score")

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

df_long_filtered <- df_long %>%
  filter(Metric != "Mean substitution standard deviation")

df_sd_fixed <- df_sd %>%
  mutate(Metric = "Mean substitution score")

maize_labels <- mapping_df_2$maize_position
names(maize_labels) <- mapping_df_2$msa_position
df_long_filtered <- df_long_filtered %>%
  left_join(mapping_df_2,
            by = c("Position" = "msa_position"))


bhlh105 <- ggplot(
  df_long,
  aes(x = Position, y = Value, fill = Metric)
) +
  
  ## --- DOMAIN REGIONS (background) ---
  geom_rect(
    data = domain_ranges_msa_2,
    aes(
      xmin = msa_start - 0.5,
      xmax = msa_end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = Domain
    ),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  
  ## --- MAIN BARS (THICKER) ---
  geom_col(
    width = 1,
    alpha = 0.6,
    position = "identity"
  ) +
  ## maize positions
  scale_x_continuous(
    breaks = mapping_df_2$msa_position,
    labels = maize_labels
  ) +
  ## --- HIGHLIGHTED BARS ---
  geom_point(
    data = df_long_filtered %>% filter(Highlight),
    aes(x = Position, y = Value),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = df_long_filtered %>% filter(Highlight),
    aes(label = maize_position),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  ## --- COLORS ---
  scale_fill_manual(values = c(
    "Mean substitution frequency" = "steelblue",
    "Mean substitution score"     = "darkred",
    "Helix-loop-helix DNA-binding domain" = "red",
    "consensus disorder prediction" = "orange"
  )) +
  
  ## --- AXES ---
  coord_cartesian(ylim = c(-2, 1), expand = FALSE) +
  
  ## --- LABELS & THEME ---
  theme_minimal(base_size = 14) +
  labs(
    title = "BHLH33 & BHLH105",
    x = "Alignment Position",
    y = "Value",
    fill = "Metric / Domain"
  ) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
bhlh105

#3 cadftr
df=read.csv("MSA_positions_table_GRMZM2G14771.csv", header = T)
# Compute consensus across each row
df_consensus <- df%>%
  rowwise() %>%
  mutate(
    Consensus = {
      # Extract all residues in this row except Position column
      residues <- c_across(-Position)
      
      # Remove gaps or missing characters
      residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # If all sequences have gap, return "-"
      if (length(residues) == 0) {
        "-"
      } else {
        # Find most frequent residue
        tbl <- table(residues)
        names(tbl)[which.max(tbl)]
      }
    }
  ) %>%
  ungroup()

# View result
head(df_consensus)
colnames(df_consensus)

# C4 species columns
C4_cols <- c("ELECO.r07.8BG0657100.1","ELECO.r07.8AG0628650.1","GRMZM2G146286_P01",
             "GRMZM2G147712","Pahal.6G178700.1","Pavag07G115100.1",
             "Urofu.6G116700.1","Sevir.6G137200.1","ELECO.r07.8BG0657100.1",
             "ELECO.r07.8AG0628650.1","GRMZM2G146286_P01","GRMZM2G147712",
             "Pahal.6G178700.1","Pavag07G115100.1","Urofu.6G116700.1","Sevir.6G137200.1")

# C3 species columns
C3_cols <- c("Bradi3g34930.1","LOC_Os08g29500.2","Chala.05G099200.1",
             "OEL18726.1","Bradi3g34930.1","LOC_Os08g29500.2",
             "Chala.05G099200.1","OEL18726.1")     


# Dataset 1: Position, Consensus, and C4 species
df_C4 <- df_consensus %>%
  select(Position, Consensus, all_of(C4_cols))

# Dataset 2: Position, Consensus, and C3 species
df_C3 <- df_consensus %>%
  select(Position, Consensus, all_of(C3_cols))


df_C3_clean <- df_C3 %>%
  filter(if_all(all_of(C3_cols), ~ !is.na(.) & . != "-"))

C3_conserved_positions <- df_C3_clean %>%
  # keep only rows where all C3 columns are identical
  filter(apply(select(., all_of(C3_cols)), 1, function(x) length(unique(x)) == 1)) %>%
  pull(Position)

df_C4_conserved_from_C3 <- df_C4 %>%
  filter(Position %in% C3_conserved_positions)

aa_cols <- setdiff(colnames(df_C4_conserved_from_C3), c("Position", "Consensus"))

n_species <- length(aa_cols)

df_consensus_mean <- df_C4_conserved_from_C3 %>%
  rowwise() %>%
  mutate(
    # collect residues for this row
    residues = list(c_across(all_of(aa_cols))),
    
    # count consensus matches
    Consensus_Matches = sum(unlist(residues) == Consensus, na.rm = TRUE),
    
    # count gaps ("-" or "")
    Gaps = sum(unlist(residues) %in% c("-", ""), na.rm = TRUE),
    
    # compute non-gap sequences
    NonGap_Count = n_species - Gaps,
    
    # compute mean consensus score
    Consensus_Mean = ifelse(
      NonGap_Count > 0,
      Consensus_Matches / NonGap_Count,
      NA_real_
    )
  ) %>%
  ungroup()  

df_consensus_mean <- df_consensus_mean %>%
  mutate(
    LowConsensus = Consensus_Mean < 0.4
  )
library(ggplot2)

ggplot(df_consensus_mean, aes(x = Position, y = Consensus_Mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  
  # Normal points
  geom_point(data = subset(df_consensus_mean, !LowConsensus),
             color = "darkblue", size = 2) +
  
  # Highlighted low-consensus positions
  geom_point(data = subset(df_consensus_mean, LowConsensus),
             color = "red", size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Consensus Mean Across Alignment (positions < 0.4 highlighted)",
    x = "Alignment Position",
    y = "Mean Consensus"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowSD)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Standard Deviation of Substitution Scores per Alignment Position",
    x = "Alignment Position",
    y = "Row SD (BLOSUM-based)"
  )
aa_cols <- C4_cols
df_row_scores <- df_consensus_mean %>%
  rowwise() %>%
  mutate(
    RowScore = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely from scoring
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Compute substitution scores only for valid AA
      scores <- sapply(valid_residues, function(aa) {
        if (aa == consensus_aa) {
          0
        } else if (aa %in% colnames(blos)) {
          blos[consensus_aa, aa]
        } else {
          NA
        }
      })
      
      # SUM the scores
      sum(scores, na.rm = TRUE)
    },
    
    RowMean = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Number of valid residues (denominator for mean)
      n_valid <- length(valid_residues)
      
      if (n_valid == 0) {
        NA   # no valid AAs at this position
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        sum(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowMean)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "black", size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution BLOSUM62 Score Across Alignment",
    x = "Alignment Position",
    y = "Row BLOSUM62 Substitution Score"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean < 0.4,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_row_scores$Consensus_Mean=1-df_row_scores$Consensus_Mean

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean > 0.75,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )

df_tmp <- df_row_scores %>%
  mutate(
    Low_Frequency = Consensus_Mean > 0.75,
    Low_Score     = RowMean < -0.5
  )
names(df_tmp)
df_tmp2 <- df_tmp
colnames(df_tmp2)[15]=c("Mean substitution frequency")
colnames(df_tmp2)[18]=c("Mean substitution score")
colnames(df_tmp2)[19]=c("Mean substitution standard deviation")

df_long <- df_tmp2 %>%
  select(
    Position,
    `Mean substitution frequency`,
    `Mean substitution score`,
    `Mean substitution standard deviation`,
    Low_Frequency,
    Low_Score
  ) %>%
  pivot_longer(
    cols = c(
      `Mean substitution frequency`,
      `Mean substitution score`,
      `Mean substitution standard deviation`
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = case_when(
      Metric == "Mean substitution frequency"          ~ Low_Frequency,
      Metric == "Mean substitution score"              ~ Low_Score,
      Metric == "Mean substitution standard deviation" ~ FALSE,
      TRUE ~ FALSE
    )
  )

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")


msa_seq <- "MDPMDIVGKSKEDVSLPKSTMVKIIKEMLPPDVRVARDAQDLLVECCV----------EFINLLSSESNEVCSREEKKTIAPEHVIKALSDLGFREYIEEVYAAYEQHKLETLDSPKAGKFTRIEMTEEEAVAEQQRMFAEARARMNNGAPKPKEPEQEPPQ-LPQ--------------------AQPQLQLHTEPQQPMQSQVQLHSQT--QH---YLQPQLQLHHQPQ--QL--PQ--VQLHSQPQL-----------Q---PQVHLHPQPQLPPQLQVHQQLQQ--P--PQVQVH------------QQPEVQPQEAQL----QSSAQQTS-QP--QPQAQLQSQ-GHSQAQLQAGL---LGQLQTQAQTG--------PDMD-S*"
maize_seq <- "MDPMDIVGKSKEDVSLPKSTMVKIIKEMLPPDVRVARDAQDLLVECCVEFINLLSSESNEVCSREEKKTIAPEHVIKALSDLGFREYIEEVYAAYEQHKLETLDSPKAGKFTRIEMTEEEAVAEQQRMFAEARARMNNGAPKPKEPEQEPPQLPQAQPQLQLHTEPQQPMQSQVQLHSQTQHYLQPQLQLHHQPQQLPQVQLHSQPQLQPQVHLHPQPQLPPQLQVHQQLQQPPQVQVHQQPEVQPQEAQLQSSAQQTSQPQPQAQLQSQGHSQAQLQAGLLGQLQTQAQTGPDMDS*"
mapping_df_2 <- map_msa_to_maize(msa_seq, maize_seq)

domains <- tribble(
  ~Domain, ~Start, ~End,
  "Histone, subunit A",	4,	143,
  "consensus disorder 2",	192,	291,
  "consensus disorder 1",	164,	183,
  #  "Zinc finger SWIM-type profile",	590,	626,
)



domain_ranges_msa_2 <- domains %>%
  rowwise() %>%
  mutate(
    msa_start = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == Start)[1] ],
    msa_end   = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == End)[1] ]
  ) %>%
  ungroup()


df_score <- df_long %>% 
  filter(Metric == "Mean substitution score")

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

df_long_filtered <- df_long %>%
  filter(Metric != "Mean substitution standard deviation")

df_sd_fixed <- df_sd %>%
  mutate(Metric = "Mean substitution score")

maize_labels <- mapping_df_2$maize_position
names(maize_labels) <- mapping_df_2$msa_position
df_long_filtered <- df_long_filtered %>%
  left_join(mapping_df_2,
            by = c("Position" = "msa_position"))


cadftr3 <- ggplot(
  df_long,
  aes(x = Position, y = Value, fill = Metric)
) +
  
  ## --- DOMAIN REGIONS (background) ---
  geom_rect(
    data = domain_ranges_msa_2,
    aes(
      xmin = msa_start - 0.5,
      xmax = msa_end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = Domain
    ),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  
  ## --- MAIN BARS (THICKER) ---
  geom_col(
    width = 1,
    alpha = 0.6,
    position = "identity"
  ) +
  ## maize positions
  scale_x_continuous(
    breaks = mapping_df_2$msa_position,
    labels = maize_labels
  ) +
  ## --- HIGHLIGHTED BARS ---
  geom_point(
    data = df_long_filtered %>% filter(Highlight),
    aes(x = Position, y = Value),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = df_long_filtered %>% filter(Highlight),
    aes(label = maize_position),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  ## --- COLORS ---
  scale_fill_manual(values = c(
    "Mean substitution frequency" = "steelblue",
    "Mean substitution score"     = "darkred",    
    "Histone, subunit A" = "red",
    "consensus disorder 1" = "orange",
    "consensus disorder 2" = "lightblue"
  )
  ) +
  
  ## --- AXES ---
  #coord_cartesian(ylim = c(-2, 2), expand = FALSE) +
  
  ## --- LABELS & THEME ---
  theme_minimal(base_size = 14) +
  labs(
    title = "CADFTR3",
    x = "Alignment Position",
    y = "Value",
    fill = "Metric / Domain"
  ) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
cadftr3

#thx8
df=read.csv("MSA_positions_table_GRMZM2G379179.csv", header = T)
# Compute consensus across each row
df_consensus <- df%>%
  rowwise() %>%
  mutate(
    Consensus = {
      # Extract all residues in this row except Position column
      residues <- c_across(-Position)
      
      # Remove gaps or missing characters
      residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # If all sequences have gap, return "-"
      if (length(residues) == 0) {
        "-"
      } else {
        # Find most frequent residue
        tbl <- table(residues)
        names(tbl)[which.max(tbl)]
      }
    }
  ) %>%
  ungroup()

# View result
head(df_consensus)
colnames(df_consensus)

# C4 species columns
C4_cols <- c("Urofu.6G176900.1","Pahal.6G242800.1","ELECO.r07.8BG0663940.1",
             "ELECO.r07.8AG0635010.1","Sevir.6G195300.1","Pavag07G170700.1","GRMZM2G379179_P01","GRMZM5G850092")

# C3 species columns
C3_cols <- c("Bradi3g38682.2","LOC_Os08g37810.1","OEL13996.1","Chala.05G054900.1")     


# Dataset 1: Position, Consensus, and C4 species
df_C4 <- df_consensus %>%
  select(Position, Consensus, all_of(C4_cols))

# Dataset 2: Position, Consensus, and C3 species
df_C3 <- df_consensus %>%
  select(Position, Consensus, all_of(C3_cols))


df_C3_clean <- df_C3 %>%
  filter(if_all(all_of(C3_cols), ~ !is.na(.) & . != "-"))

C3_conserved_positions <- df_C3_clean %>%
  # keep only rows where all C3 columns are identical
  filter(apply(select(., all_of(C3_cols)), 1, function(x) length(unique(x)) == 1)) %>%
  pull(Position)

df_C4_conserved_from_C3 <- df_C4 %>%
  filter(Position %in% C3_conserved_positions)

aa_cols <- setdiff(colnames(df_C4_conserved_from_C3), c("Position", "Consensus"))

n_species <- length(aa_cols)

df_consensus_mean <- df_C4_conserved_from_C3 %>%
  rowwise() %>%
  mutate(
    # collect residues for this row
    residues = list(c_across(all_of(aa_cols))),
    
    # count consensus matches
    Consensus_Matches = sum(unlist(residues) == Consensus, na.rm = TRUE),
    
    # count gaps ("-" or "")
    Gaps = sum(unlist(residues) %in% c("-", ""), na.rm = TRUE),
    
    # compute non-gap sequences
    NonGap_Count = n_species - Gaps,
    
    # compute mean consensus score
    Consensus_Mean = ifelse(
      NonGap_Count > 0,
      Consensus_Matches / NonGap_Count,
      NA_real_
    )
  ) %>%
  ungroup()  

df_consensus_mean <- df_consensus_mean %>%
  mutate(
    LowConsensus = Consensus_Mean < 0.4
  )
library(ggplot2)

ggplot(df_consensus_mean, aes(x = Position, y = Consensus_Mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  
  # Normal points
  geom_point(data = subset(df_consensus_mean, !LowConsensus),
             color = "darkblue", size = 2) +
  
  # Highlighted low-consensus positions
  geom_point(data = subset(df_consensus_mean, LowConsensus),
             color = "red", size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Consensus Mean Across Alignment (positions < 0.4 highlighted)",
    x = "Alignment Position",
    y = "Mean Consensus"
  )

ggplot(df_row_scores, aes(x = Position, y = RowSD)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Standard Deviation of Substitution Scores per Alignment Position",
    x = "Alignment Position",
    y = "Row SD (BLOSUM-based)"
  )
aa_cols <- C4_cols
df_row_scores <- df_consensus_mean %>%
  rowwise() %>%
  mutate(
    RowScore = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely from scoring
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Compute substitution scores only for valid AA
      scores <- sapply(valid_residues, function(aa) {
        if (aa == consensus_aa) {
          0
        } else if (aa %in% colnames(blos)) {
          blos[consensus_aa, aa]
        } else {
          NA
        }
      })
      
      # SUM the scores
      sum(scores, na.rm = TRUE)
    },
    
    RowMean = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Number of valid residues (denominator for mean)
      n_valid <- length(valid_residues)
      
      if (n_valid == 0) {
        NA   # no valid AAs at this position
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        sum(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowMean)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "black", size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution BLOSUM62 Score Across Alignment",
    x = "Alignment Position",
    y = "Row BLOSUM62 Substitution Score"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean < 0.4,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_row_scores$Consensus_Mean=1-df_row_scores$Consensus_Mean

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean > 0.75,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_tmp <- df_row_scores %>%
  mutate(
    Low_Frequency = Consensus_Mean > 0.75,
    Low_Score     = RowMean < -0.5
  )
names(df_tmp)
df_tmp2 <- df_tmp
colnames(df_tmp2)[15]=c("Mean substitution frequency")
colnames(df_tmp2)[18]=c("Mean substitution score")
colnames(df_tmp2)[19]=c("Mean substitution standard deviation")

df_long <- df_tmp2 %>%
  select(
    Position,
    `Mean substitution frequency`,
    `Mean substitution score`,
    `Mean substitution standard deviation`,
    Low_Frequency,
    Low_Score
  ) %>%
  pivot_longer(
    cols = c(
      `Mean substitution frequency`,
      `Mean substitution score`,
      `Mean substitution standard deviation`
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = case_when(
      Metric == "Mean substitution frequency"          ~ Low_Frequency,
      Metric == "Mean substitution score"              ~ Low_Score,
      Metric == "Mean substitution standard deviation" ~ FALSE,
      TRUE ~ FALSE
    )
  )

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")


msa_seq <- "M--------A-EEP---PPP-----------PL-ARKSAPAQPWSHVETTHLIDAYEERWTALRRGQLKAHQWEEVAAEVAARCAATPGVVAQRKTGTQCRHKLEKLRKRYRTEGARPVTSLWPYFRRMDRLERGPLAVASSAYPA--ATGSPPAADGD--EEGEE-EEEEVEEEEEEEDVQEENE-------EEEELAPRNNNTRSINGIIREFG------TG---LAPRHPQLQLHP-----PPPSITPSIAPPRKRVAYEAFQAKAAAAAV---KAKDDEEEA--MEM-AR-RRGSSGRPVAQLSAVLRDFGEGVMRLERRRMEVQWEIERGWQEADARHARMLQDAQRQLRDTVAG-A-CALPPKKARRDHG---DS-*-------------------------------------"
maize_seq <- "MAEEPPPPPLARKSAPAQPWSHVETTHLIDAYEERWTALRRGQLKAHQWEEVAAEVAARCAATPGVVAQRKTGTQCRHKLEKLRKRYRTEGARPVTSLWPYFRRMDRLERGPLAVASSAYPAATGSPPAADGDEEGEEEEEEVEEEEEEEDVQEENEEEEELAPRNNNTRSINGIIREFGTGLAPRHPQLQLHPPPPSITPSIAPPRKRVAYEAFQAKAAAAAVKAKDDEEEAMEMARRRGSSGRPVAQLSAVLRDFGEGVMRLERRRMEVQWEIERGWQEADARHARMLQDAQRQLRDTVAGACALPPKKARRDHGDS*"
mapping_df_2 <- map_msa_to_maize(msa_seq, maize_seq)

domains <- tribble(
  ~Domain, ~Start, ~End,
  "Myb/SANT-like DNA-binding domain",	19,	108,
  "consensus disorder prediction",	116,	205,
  #  "Zinc finger SWIM-type profile",	590,	626,
)




domain_ranges_msa_2 <- domains %>%
  rowwise() %>%
  mutate(
    msa_start = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == Start)[1] ],
    msa_end   = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == End)[1] ]
  ) %>%
  ungroup()


df_score <- df_long %>% 
  filter(Metric == "Mean substitution score")

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

df_long_filtered <- df_long %>%
  filter(Metric != "Mean substitution standard deviation")

df_sd_fixed <- df_sd %>%
  mutate(Metric = "Mean substitution score")

maize_labels <- mapping_df_2$maize_position
names(maize_labels) <- mapping_df_2$msa_position
df_long_filtered <- df_long_filtered %>%
  left_join(mapping_df_2,
            by = c("Position" = "msa_position"))


thx8 <- ggplot(
  df_long,
  aes(x = Position, y = Value, fill = Metric)
) +
  
  ## --- DOMAIN REGIONS (background) ---
  geom_rect(
    data = domain_ranges_msa_2,
    aes(
      xmin = msa_start - 0.5,
      xmax = msa_end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = Domain
    ),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  
  ## --- MAIN BARS (THICKER) ---
  geom_col(
    width = 1,
    alpha = 0.6,
    position = "identity"
  ) +
  ## maize positions
  scale_x_continuous(
    breaks = mapping_df_2$msa_position,
    labels = maize_labels
  ) +
  ## --- HIGHLIGHTED BARS ---
  geom_point(
    data = df_long_filtered %>% filter(Highlight),
    aes(x = Position, y = Value),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = df_long_filtered %>% filter(Highlight),
    aes(label = maize_position),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  ## --- COLORS ---
  scale_fill_manual(values = c(
    "Mean substitution frequency" = "steelblue",
    "Mean substitution score"     = "darkred",
    "Myb/SANT-like DNA-binding domain" = "red",
    "consensus disorder prediction" = "orange"
  )) +
  ## --- AXES ---
  #coord_cartesian(ylim = c(-2, 2), expand = FALSE) +
  
  ## --- LABELS & THEME ---
  theme_minimal(base_size = 14) +
  labs(
    title = "THX8",
    x = "Alignment Position",
    y = "Value",
    fill = "Metric / Domain"
  ) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
thx8
#GRMZM2G005155
setwd("~/angeo_C4/bootstrap/GRMZM2G005155")
df=read.csv("MSA_positions_table.csv", header = T)
# Compute consensus across each row
df_consensus <- df%>%
  rowwise() %>%
  mutate(
    Consensus = {
      # Extract all residues in this row except Position column
      residues <- c_across(-Position)
      
      # Remove gaps or missing characters
      residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # If all sequences have gap, return "-"
      if (length(residues) == 0) {
        "-"
      } else {
        # Find most frequent residue
        tbl <- table(residues)
        names(tbl)[which.max(tbl)]
      }
    }
  ) %>%
  ungroup()

# View result
head(df_consensus)
colnames(df_consensus)

# C4 species columns
C4_cols <- c("GRMZM2G148220","GRMZM2G128953_P01","GRMZM2G137387_P01","ELECO.r07.4BG0334990.1","ELECO.r07.4AG0304100.1",
             "GRMZM2G137289_P01","Sevir.7G247100.1","Pahal.7G290300.1","Urofu.7G275200.1","Pavag06G226400.1",
             "GRMZM2G005155_P03","GRMZM2G135018_P01","GRMZM2G072997_P01")

# C3 species columns
C3_cols <- c("OEL20795.1","Chala.06G231700.1","Bradi5g21700.2","LOC_Os04g52410.2")     


# Dataset 1: Position, Consensus, and C4 species
df_C4 <- df_consensus %>%
  select(Position, Consensus, all_of(C4_cols))

# Dataset 2: Position, Consensus, and C3 species
df_C3 <- df_consensus %>%
  select(Position, Consensus, all_of(C3_cols))


df_C3_clean <- df_C3 %>%
  filter(if_all(all_of(C3_cols), ~ !is.na(.) & . != "-"))

C3_conserved_positions <- df_C3_clean %>%
  # keep only rows where all C3 columns are identical
  filter(apply(select(., all_of(C3_cols)), 1, function(x) length(unique(x)) == 1)) %>%
  pull(Position)

df_C4_conserved_from_C3 <- df_C4 %>%
  filter(Position %in% C3_conserved_positions)

aa_cols <- setdiff(colnames(df_C4_conserved_from_C3), c("Position", "Consensus"))

n_species <- length(aa_cols)

df_consensus_mean <- df_C4_conserved_from_C3 %>%
  rowwise() %>%
  mutate(
    # collect residues for this row
    residues = list(c_across(all_of(aa_cols))),
    
    # count consensus matches
    Consensus_Matches = sum(unlist(residues) == Consensus, na.rm = TRUE),
    
    # count gaps ("-" or "")
    Gaps = sum(unlist(residues) %in% c("-", ""), na.rm = TRUE),
    
    # compute non-gap sequences
    NonGap_Count = n_species - Gaps,
    
    # compute mean consensus score
    Consensus_Mean = ifelse(
      NonGap_Count > 0,
      Consensus_Matches / NonGap_Count,
      NA_real_
    )
  ) %>%
  ungroup()  

df_consensus_mean <- df_consensus_mean %>%
  mutate(
    LowConsensus = Consensus_Mean < 0.4
  )
library(ggplot2)

ggplot(df_consensus_mean, aes(x = Position, y = Consensus_Mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  
  # Normal points
  geom_point(data = subset(df_consensus_mean, !LowConsensus),
             color = "darkblue", size = 2) +
  
  # Highlighted low-consensus positions
  geom_point(data = subset(df_consensus_mean, LowConsensus),
             color = "red", size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Consensus Mean Across Alignment (positions < 0.4 highlighted)",
    x = "Alignment Position",
    y = "Mean Consensus"
  )

ggplot(df_row_scores, aes(x = Position, y = RowSD)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Standard Deviation of Substitution Scores per Alignment Position",
    x = "Alignment Position",
    y = "Row SD (BLOSUM-based)"
  )
aa_cols <- C4_cols
df_row_scores <- df_consensus_mean %>%
  rowwise() %>%
  mutate(
    RowScore = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely from scoring
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Compute substitution scores only for valid AA
      scores <- sapply(valid_residues, function(aa) {
        if (aa == consensus_aa) {
          0
        } else if (aa %in% colnames(blos)) {
          blos[consensus_aa, aa]
        } else {
          NA
        }
      })
      
      # SUM the scores
      sum(scores, na.rm = TRUE)
    },
    
    RowMean = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Number of valid residues (denominator for mean)
      n_valid <- length(valid_residues)
      
      if (n_valid == 0) {
        NA   # no valid AAs at this position
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        sum(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowMean)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "black", size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution BLOSUM62 Score Across Alignment",
    x = "Alignment Position",
    y = "Row BLOSUM62 Substitution Score"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean < 0.4,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_row_scores$Consensus_Mean=1-df_row_scores$Consensus_Mean

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean > 0.75,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )

df_tmp <- df_row_scores %>%
  mutate(
    Low_Frequency = Consensus_Mean > 0.75,
    Low_Score     = RowMean < -0.5
  )
names(df_tmp)
df_tmp2 <- df_tmp
colnames(df_tmp2)[20]=c("Mean substitution frequency")
colnames(df_tmp2)[23]=c("Mean substitution score")
colnames(df_tmp2)[24]=c("Mean substitution standard deviation")

df_long <- df_tmp2 %>%
  select(
    Position,
    `Mean substitution frequency`,
    `Mean substitution score`,
    `Mean substitution standard deviation`,
    Low_Frequency,
    Low_Score
  ) %>%
  pivot_longer(
    cols = c(
      `Mean substitution frequency`,
      `Mean substitution score`,
      `Mean substitution standard deviation`
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = case_when(
      Metric == "Mean substitution frequency"          ~ Low_Frequency,
      Metric == "Mean substitution score"              ~ Low_Score,
      Metric == "Mean substitution standard deviation" ~ FALSE,
      TRUE ~ FALSE
    )
  )

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

msa_seq <- "MGRGKVELKKIENPTNRQVTFSKRRMGLFKKANELAILCDAQIGVIIFSGSGRMYEYSSPPWRIASVFDRYLKAPSTRFEEMDIQQKIVQEMTRMKDERNRLRMIMAQYMAEDLASFSAQDLSNLEQQIEFSLYKVRLRKQELLDQQLLEIHQREMHMPAEQGGYLCLMNPASGQHQQAGEMVNPRPFPWWDVGASGSLLHGRDAESSMTALGLSPQLHGYRLQPRQPNLQDADIHGWL"
maize_seq <- "MGRGKVELKKIENPTNRQVTFSKRRMGLFKKANELAILCDAQIGVIIFSGSGRMYEYSSPPWRIASVFDRYLKAPSTRFEEMDIQQKIVQEMTRMKDERNRLRMIMAQYMAEDLASFSAQDLSNLEQQIEFSLYKVRLRKQELLDQQLLEIHQREMHMPAEQGGYLCLMNPAAAIASGQHQQAGEMVGINPRPFPWWDVGASGSGSGSQSQQQLLHGRDAAESSMTALGLSPQLHGYRLQPRQPNLQQDADIHGWL*"
mapping_df_2 <- map_msa_to_maize(msa_seq, maize_seq)

domains <- tribble(
  ~Domain, ~Start, ~End,
  "SRF-like domain",	2,	79,
  "K_BOX domain",	85,	182,
  #  "Zinc finger SWIM-type profile",	590,	626,
)


domain_ranges_msa_2 <- domains %>%
  rowwise() %>%
  mutate(
    msa_start = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == Start)[1] ],
    msa_end   = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == End)[1] ]
  ) %>%
  ungroup()
domain_ranges_msa_2

df_score <- df_long %>% 
  filter(Metric == "Mean substitution score")

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

df_long_filtered <- df_long %>%
  filter(Metric != "Mean substitution standard deviation")

df_sd_fixed <- df_sd %>%
  mutate(Metric = "Mean substitution score")

maize_labels <- mapping_df_2$maize_position
names(maize_labels) <- mapping_df_2$msa_position
df_long_filtered <- df_long_filtered %>%
  left_join(mapping_df_2,
            by = c("Position" = "msa_position"))


mads9 <- ggplot(
  df_long,
  aes(x = Position, y = Value, fill = Metric)
) +
  
  ## --- DOMAIN REGIONS (background) ---
  geom_rect(
    data = domain_ranges_msa_2,
    aes(
      xmin = msa_start - 0.5,
      xmax = msa_end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = Domain
    ),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  
  ## --- MAIN BARS (THICKER) ---
  geom_col(
    width = 1,
    alpha = 0.6,
    position = "identity"
  ) +
  ## maize positions
  scale_x_continuous(
    breaks = mapping_df_2$msa_position,
    labels = maize_labels
  ) +
  ## --- HIGHLIGHTED BARS ---
  geom_point(
    data = df_long_filtered %>% filter(Highlight),
    aes(x = Position, y = Value),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = df_long_filtered %>% filter(Highlight),
    aes(label = maize_position),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  ## --- COLORS ---
  scale_fill_manual(
    values = c(
      "Mean substitution frequency" = "steelblue",
      "Mean substitution score"     = "darkred",
      "SRF-like domain"     = "red",
      "K_BOX domain" = "orange"
    )
  ) +
  
  ## --- AXES ---
  #coord_cartesian(ylim = c(-2, 2), expand = FALSE) +
  
  ## --- LABELS & THEME ---
  theme_minimal(base_size = 14) +
  labs(
    title = "MADS9",
    x = "Alignment Position",
    y = "Value",
    fill = "Metric / Domain"
  ) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
mads9

#GRMZM2G036837
setwd("~/angeo_C4/bootstrap/GRMZM2G036837")
df=read.csv("MSA_positions_table.csv", header = T)
# Compute consensus across each row
df_consensus <- df%>%
  rowwise() %>%
  mutate(
    Consensus = {
      # Extract all residues in this row except Position column
      residues <- c_across(-Position)
      
      # Remove gaps or missing characters
      residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # If all sequences have gap, return "-"
      if (length(residues) == 0) {
        "-"
      } else {
        # Find most frequent residue
        tbl <- table(residues)
        names(tbl)[which.max(tbl)]
      }
    }
  ) %>%
  ungroup()

# View result
head(df_consensus)
colnames(df_consensus)

# C4 species columns
C4_cols <- c("GRMZM5G871727","Sevir.5G055700.1","ELECO.r07.1AG0005070.1","ELECO.r07.1BG0053580.1",
             "GRMZM2G036837_P02","GRMZM2G086994_P02","Pahal.5G426900.1","Urofu.5G383800.1","Pavag03G106600.1")

# C3 species columns
C3_cols <- c("Bradi2g09350.1","Chala.08G157200.1","LOC_Os01g15460.1")     


# Dataset 1: Position, Consensus, and C4 species
df_C4 <- df_consensus %>%
  select(Position, Consensus, all_of(C4_cols))

# Dataset 2: Position, Consensus, and C3 species
df_C3 <- df_consensus %>%
  select(Position, Consensus, all_of(C3_cols))


df_C3_clean <- df_C3 %>%
  filter(if_all(all_of(C3_cols), ~ !is.na(.) & . != "-"))

C3_conserved_positions <- df_C3_clean %>%
  # keep only rows where all C3 columns are identical
  filter(apply(select(., all_of(C3_cols)), 1, function(x) length(unique(x)) == 1)) %>%
  pull(Position)

df_C4_conserved_from_C3 <- df_C4 %>%
  filter(Position %in% C3_conserved_positions)

aa_cols <- setdiff(colnames(df_C4_conserved_from_C3), c("Position", "Consensus"))

n_species <- length(aa_cols)

df_consensus_mean <- df_C4_conserved_from_C3 %>%
  rowwise() %>%
  mutate(
    # collect residues for this row
    residues = list(c_across(all_of(aa_cols))),
    
    # count consensus matches
    Consensus_Matches = sum(unlist(residues) == Consensus, na.rm = TRUE),
    
    # count gaps ("-" or "")
    Gaps = sum(unlist(residues) %in% c("-", ""), na.rm = TRUE),
    
    # compute non-gap sequences
    NonGap_Count = n_species - Gaps,
    
    # compute mean consensus score
    Consensus_Mean = ifelse(
      NonGap_Count > 0,
      Consensus_Matches / NonGap_Count,
      NA_real_
    )
  ) %>%
  ungroup()  

df_consensus_mean <- df_consensus_mean %>%
  mutate(
    LowConsensus = Consensus_Mean < 0.4
  )
library(ggplot2)

ggplot(df_consensus_mean, aes(x = Position, y = Consensus_Mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  
  # Normal points
  geom_point(data = subset(df_consensus_mean, !LowConsensus),
             color = "darkblue", size = 2) +
  
  # Highlighted low-consensus positions
  geom_point(data = subset(df_consensus_mean, LowConsensus),
             color = "red", size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Consensus Mean Across Alignment (positions < 0.4 highlighted)",
    x = "Alignment Position",
    y = "Mean Consensus"
  )

ggplot(df_row_scores, aes(x = Position, y = RowSD)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Standard Deviation of Substitution Scores per Alignment Position",
    x = "Alignment Position",
    y = "Row SD (BLOSUM-based)"
  )
aa_cols <- C4_cols
df_row_scores <- df_consensus_mean %>%
  rowwise() %>%
  mutate(
    RowScore = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely from scoring
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Compute substitution scores only for valid AA
      scores <- sapply(valid_residues, function(aa) {
        if (aa == consensus_aa) {
          0
        } else if (aa %in% colnames(blos)) {
          blos[consensus_aa, aa]
        } else {
          NA
        }
      })
      
      # SUM the scores
      sum(scores, na.rm = TRUE)
    },
    
    RowMean = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Number of valid residues (denominator for mean)
      n_valid <- length(valid_residues)
      
      if (n_valid == 0) {
        NA   # no valid AAs at this position
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        sum(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowMean)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "black", size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution BLOSUM62 Score Across Alignment",
    x = "Alignment Position",
    y = "Row BLOSUM62 Substitution Score"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean < 0.4,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_row_scores$Consensus_Mean=1-df_row_scores$Consensus_Mean

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean > 0.75,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )

df_tmp <- df_row_scores %>%
  mutate(
    Low_Frequency = Consensus_Mean > 0.75,
    Low_Score     = RowMean < -0.5
  )
names(df_tmp)
df_tmp2 <- df_tmp
colnames(df_tmp2)[16]=c("Mean substitution frequency")
colnames(df_tmp2)[19]=c("Mean substitution score")
colnames(df_tmp2)[20]=c("Mean substitution standard deviation")

df_long <- df_tmp2 %>%
  select(
    Position,
    `Mean substitution frequency`,
    `Mean substitution score`,
    `Mean substitution standard deviation`,
    Low_Frequency,
    Low_Score
  ) %>%
  pivot_longer(
    cols = c(
      `Mean substitution frequency`,
      `Mean substitution score`,
      `Mean substitution standard deviation`
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = case_when(
      Metric == "Mean substitution frequency"          ~ Low_Frequency,
      Metric == "Mean substitution score"              ~ Low_Score,
      Metric == "Mean substitution standard deviation" ~ FALSE,
      TRUE ~ FALSE
    )
  )

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

msa_seq <- "MEPHAEQAVAAVAGGEGGTASPGTGLEGPVLRRGLDGGGEGEDGEGAQEANARLPERPGEADCGYYLRTGACGFGERCRYNHPRDRGGTEFGGGAKNGAAQDFPERQGQPVCEYYLKTGTCKFGSNCKYHHPKQDGSVQSVILNNNGFPLRLGEKECSYYMKTGQCKFGSTCKFHHPEFGGIPVTPGIYPPLQSPSVPSPHTYAP--NWQMGRSPAVPGSYIPGSYTPMMISSGMVPLQGWSPYPASVNPVASGGAQQTVQAGPLYGIGHHGSSTAIAYGGTYLPYSSSAGQSSNNHQEHGFPERPGQPECQYYMRTGDCKFGTTCKYNHPRDWSTPKSNYMFSHLCLPLRPGAQPCAYYAQNGYCRYGVACKYDHSMGTLGYSSSALPLSDMPIAPYPISFSVATLAPSSSSPEYISTKDPSINHVVSPVAGPAPVGAILPKGVFHHDTIMQTQTPTS-AGSSSPGGGR"
maize_seq <- "MEPPHAEQAVAAVAAGGEGGTASPGTGLEGPVLRRGLDGGGEGEDGELGAGQEANARLPERPGEADCGYYLRTGACGFGERCRYNHPRDRGGTEFGGGAKNGAAQDFPERQGQPVCEYYLKTGTCKFGSNCKYHHPKQDGSVQSVILNNNGFPLRLGEKECSYYMKTGQCKFGSTCKFHHPEFGGIPVTPGIYPPLQSPSVPSPHTYAPNWQMGRSPAVPGSYIPGSYTPMMISSGMVPLQGWSPYPASVNPVASGGAQQTVQAGPLYGIGHHGSSTAIAYGGTYLPYSSSAGQSSNNHQEHGFPERPGQPECQYYMRTGDCKFGTTCKYNHPRDWSTPKSNYMFSHLCLPLRPGAQPCAYYAQNGYCRYGVACKYDHSMGTLGYSSSALPLSDMPIAPYPISFSVATLAPSSSSPEYISTKDPSINHVVSPVAGPAPVGAILPKGVFHHDTIMQTQTPTSAGSSSPGGGR*"
mapping_df_2 <- map_msa_to_maize(msa_seq, maize_seq)

domains <- tribble(
  ~Domain, ~Start, ~End,
  "Zinc finger C3H1-type profile 1",	349,	377,
  "Zinc finger C3H1-type profile 2",	106,	134,
  "CCCH zinc finger 1",	307,	331,
  "CCCH zinc finger 2",	153,	178,
  #  "Zinc finger SWIM-type profile",	590,	626,
)


domain_ranges_msa_2 <- domains %>%
  rowwise() %>%
  mutate(
    msa_start = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == Start)[1] ],
    msa_end   = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == End)[1] ]
  ) %>%
  ungroup()


df_score <- df_long %>% 
  filter(Metric == "Mean substitution score")

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

df_long_filtered <- df_long %>%
  filter(Metric != "Mean substitution standard deviation")

df_sd_fixed <- df_sd %>%
  mutate(Metric = "Mean substitution score")

maize_labels <- mapping_df_2$maize_position
names(maize_labels) <- mapping_df_2$msa_position
df_long_filtered <- df_long_filtered %>%
  left_join(mapping_df_2,
            by = c("Position" = "msa_position"))


c3h28 <- ggplot(
  df_long,
  aes(x = Position, y = Value, fill = Metric)
) +
  
  ## --- DOMAIN REGIONS (background) ---
  geom_rect(
    data = domain_ranges_msa_2,
    aes(
      xmin = msa_start - 0.5,
      xmax = msa_end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = Domain
    ),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  
  ## --- MAIN BARS (THICKER) ---
  geom_col(
    width = 1,
    alpha = 0.6,
    position = "identity"
  ) +
  ## maize positions
  scale_x_continuous(
    breaks = mapping_df_2$msa_position,
    labels = maize_labels
  ) +
  ## --- HIGHLIGHTED BARS ---
  geom_point(
    data = df_long_filtered %>% filter(Highlight),
    aes(x = Position, y = Value),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = df_long_filtered %>% filter(Highlight),
    aes(label = maize_position),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  ## --- COLORS ---
  scale_fill_manual(
    values = c(
      "Mean substitution frequency" = "steelblue",
      "Mean substitution score"     = "darkred",
      "Zinc finger C3H1-type profile 1"     = "red",
      "Zinc finger C3H1-type profile 2"= "orange",
      "CCCH zinc finger 1" = "purple",
      "CCCH zinc finger 2" = "yellow" 
    )
  ) +
  
  ## --- AXES ---
  #coord_cartesian(ylim = c(-2, 2), expand = FALSE) +
  
  ## --- LABELS & THEME ---
  theme_minimal(base_size = 14) +
  labs(
    title = "C3H28",
    x = "Alignment Position",
    y = "Value",
    fill = "Metric / Domain"
  ) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
c3h28
#GRMZM2G042895
setwd("~/angeo_C4/bootstrap/GRMZM2G042895")
df=read.csv("MSA_positions_table.csv", header = T)
# Compute consensus across each row
df_consensus <- df%>%
  rowwise() %>%
  mutate(
    Consensus = {
      # Extract all residues in this row except Position column
      residues <- c_across(-Position)
      
      # Remove gaps or missing characters
      residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # If all sequences have gap, return "-"
      if (length(residues) == 0) {
        "-"
      } else {
        # Find most frequent residue
        tbl <- table(residues)
        names(tbl)[which.max(tbl)]
      }
    }
  ) %>%
  ungroup()

# View result
head(df_consensus)
colnames(df_consensus)

# C4 species columns
C4_cols <- c("GRMZM2G042895","Pavag06G047800.1","Pahal.7G091500.1",
             "Urofu.7G060000.1","Sevir.7G050300.1.p","Sevir.7G050400.1",
             "Pahal.7G091300.1.p","ELECO.r07.4BG0352570.1","ELECO.r07.4AG0323410.1")

# C3 species columns
C3_cols <- c("Bradi5g05110.1","Chala.06G054300.1.p","LOC_Os04g23550.1","Chala.06G054400.1")     


# Dataset 1: Position, Consensus, and C4 species
df_C4 <- df_consensus %>%
  select(Position, Consensus, all_of(C4_cols))

# Dataset 2: Position, Consensus, and C3 species
df_C3 <- df_consensus %>%
  select(Position, Consensus, all_of(C3_cols))


df_C3_clean <- df_C3 %>%
  filter(if_all(all_of(C3_cols), ~ !is.na(.) & . != "-"))

C3_conserved_positions <- df_C3_clean %>%
  # keep only rows where all C3 columns are identical
  filter(apply(select(., all_of(C3_cols)), 1, function(x) length(unique(x)) == 1)) %>%
  pull(Position)

df_C4_conserved_from_C3 <- df_C4 %>%
  filter(Position %in% C3_conserved_positions)

aa_cols <- setdiff(colnames(df_C4_conserved_from_C3), c("Position", "Consensus"))

n_species <- length(aa_cols)

df_consensus_mean <- df_C4_conserved_from_C3 %>%
  rowwise() %>%
  mutate(
    # collect residues for this row
    residues = list(c_across(all_of(aa_cols))),
    
    # count consensus matches
    Consensus_Matches = sum(unlist(residues) == Consensus, na.rm = TRUE),
    
    # count gaps ("-" or "")
    Gaps = sum(unlist(residues) %in% c("-", ""), na.rm = TRUE),
    
    # compute non-gap sequences
    NonGap_Count = n_species - Gaps,
    
    # compute mean consensus score
    Consensus_Mean = ifelse(
      NonGap_Count > 0,
      Consensus_Matches / NonGap_Count,
      NA_real_
    )
  ) %>%
  ungroup()  

df_consensus_mean <- df_consensus_mean %>%
  mutate(
    LowConsensus = Consensus_Mean < 0.4
  )
library(ggplot2)

ggplot(df_consensus_mean, aes(x = Position, y = Consensus_Mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  
  # Normal points
  geom_point(data = subset(df_consensus_mean, !LowConsensus),
             color = "darkblue", size = 2) +
  
  # Highlighted low-consensus positions
  geom_point(data = subset(df_consensus_mean, LowConsensus),
             color = "red", size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Consensus Mean Across Alignment (positions < 0.4 highlighted)",
    x = "Alignment Position",
    y = "Mean Consensus"
  )

ggplot(df_row_scores, aes(x = Position, y = RowSD)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Standard Deviation of Substitution Scores per Alignment Position",
    x = "Alignment Position",
    y = "Row SD (BLOSUM-based)"
  )
aa_cols <- C4_cols
df_row_scores <- df_consensus_mean %>%
  rowwise() %>%
  mutate(
    RowScore = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely from scoring
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Compute substitution scores only for valid AA
      scores <- sapply(valid_residues, function(aa) {
        if (aa == consensus_aa) {
          0
        } else if (aa %in% colnames(blos)) {
          blos[consensus_aa, aa]
        } else {
          NA
        }
      })
      
      # SUM the scores
      sum(scores, na.rm = TRUE)
    },
    
    RowMean = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Number of valid residues (denominator for mean)
      n_valid <- length(valid_residues)
      
      if (n_valid == 0) {
        NA   # no valid AAs at this position
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        sum(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowMean)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "black", size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution BLOSUM62 Score Across Alignment",
    x = "Alignment Position",
    y = "Row BLOSUM62 Substitution Score"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean < 0.4,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_row_scores$Consensus_Mean=1-df_row_scores$Consensus_Mean

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean > 0.75,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )

df_tmp <- df_row_scores %>%
  mutate(
    Low_Frequency = Consensus_Mean > 0.75,
    Low_Score     = RowMean < -0.5
  )
names(df_tmp)
df_tmp2 <- df_tmp
colnames(df_tmp2)[16]=c("Mean substitution frequency")
colnames(df_tmp2)[19]=c("Mean substitution score")
colnames(df_tmp2)[20]=c("Mean substitution standard deviation")

df_long <- df_tmp2 %>%
  select(
    Position,
    `Mean substitution frequency`,
    `Mean substitution score`,
    `Mean substitution standard deviation`,
    Low_Frequency,
    Low_Score
  ) %>%
  pivot_longer(
    cols = c(
      `Mean substitution frequency`,
      `Mean substitution score`,
      `Mean substitution standard deviation`
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = case_when(
      Metric == "Mean substitution frequency"          ~ Low_Frequency,
      Metric == "Mean substitution score"              ~ Low_Score,
      Metric == "Mean substitution standard deviation" ~ FALSE,
      TRUE ~ FALSE
    )
  )

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

msa_seq <- "MEMGDSYYWEMQQYLESEELSLYMGTQDDALSCYDSSSPDGSISNSWAPAGVSREGGAAAANKNILMERDRRRKLNEKLYALRSVVPNITKMDKASIIKDAIEYIEQLQAEERRALQALEAGEGARCHGHGEEARVVLQQ------------PA-AAPAPVEVLELRVSEVGDRVLVVNVTCSKGRDAMARVCRAVEELRLRVITASVTSVAGCLMHTIFVEVDQTNRIQIKHMIEAALAQLDDSASPPSVMSY"
maize_seq <- "MEMGDSFEYYWEMQQYLESEELSLYMGTQDDALSCYDSSSPDGSISNSSWAPAGVAATASEKREGPGGAAAANKNILMERDRRRKLNEKLYALRSVVPNITKMDKASIIKDAIEYIEQLQAEERRALQALEAGEGARCGGHGHGEEARVVLQQPAAAPAPVEVLELRVSEVGDRVLVVNVTCSKGRDAMARVCRAVEELRLRVITASVTSVAGCLMHTIFVEVDSDQTNRIQIKHMIEAALAQLDDASASPPSVMSYY*"
mapping_df_2 <- map_msa_to_maize(msa_seq, maize_seq)

domains <- tribble(
  ~Domain, ~Start, ~End,
  "HLH, helix-loop-helix domain",	70,	120,
  "bHLH-TF_ACT-like_plant",	165,	241,
  #  "Zinc finger SWIM-type profile",	590,	626,
)


domain_ranges_msa_2 <- domains %>%
  rowwise() %>%
  mutate(
    msa_start = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == Start)[1] ],
    msa_end   = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == End)[1] ]
  ) %>%
  ungroup()


df_score <- df_long %>% 
  filter(Metric == "Mean substitution score")

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

df_long_filtered <- df_long %>%
  filter(Metric != "Mean substitution standard deviation")

df_sd_fixed <- df_sd %>%
  mutate(Metric = "Mean substitution score")

maize_labels <- mapping_df_2$maize_position
names(maize_labels) <- mapping_df_2$msa_position
df_long_filtered <- df_long_filtered %>%
  left_join(mapping_df_2,
            by = c("Position" = "msa_position"))


bhlh116 <- ggplot(
  df_long,
  aes(x = Position, y = Value, fill = Metric)
) +
  
  ## --- DOMAIN REGIONS (background) ---
  geom_rect(
    data = domain_ranges_msa_2,
    aes(
      xmin = msa_start - 0.5,
      xmax = msa_end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = Domain
    ),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  
  ## --- MAIN BARS (THICKER) ---
  geom_col(
    width = 1,
    alpha = 0.6,
    position = "identity"
  )  +
  ## maize positions
  scale_x_continuous(
    breaks = mapping_df_2$msa_position,
    labels = maize_labels
  ) +
  ## --- HIGHLIGHTED BARS ---
  geom_point(
    data = df_long_filtered %>% filter(Highlight),
    aes(x = Position, y = Value),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = df_long_filtered %>% filter(Highlight),
    aes(label = maize_position),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  ## --- COLORS ---
  scale_fill_manual(
    values = c(
      "Mean substitution frequency" = "steelblue",
      "Mean substitution score"     = "darkred",
      "HLH, helix-loop-helix domain"     = "red",
      "bHLH-TF_ACT-like_plant"= "orange" 
    )
  ) +
  
  ## --- AXES ---
  #coord_cartesian(ylim = c(-2, 2), expand = FALSE) +
  
  ## --- LABELS & THEME ---
  theme_minimal(base_size = 14) +
  labs(
    title = "BHLH116",
    x = "Alignment Position",
    y = "Value",
    fill = "Metric / Domain"
  ) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
bhlh116


#scarecrow
setwd("~/angeo_C4/positivecontrol")
df=read.csv("MSA_positions_table_GRMZM2G015080.csv", header = T)
# Compute consensus across each row
df_consensus <- df%>%
  rowwise() %>%
  mutate(
    Consensus = {
      # Extract all residues in this row except Position column
      residues <- c_across(-Position)
      
      # Remove gaps or missing characters
      residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # If all sequences have gap, return "-"
      if (length(residues) == 0) {
        "-"
      } else {
        # Find most frequent residue
        tbl <- table(residues)
        names(tbl)[which.max(tbl)]
      }
    }
  ) %>%
  ungroup()

# View result
head(df_consensus)
colnames(df_consensus)

# C4 species columns
C4_cols <- c("ELECO.r07.9BG0696540.1","ELECO.r07.9AG0672860.1","Urofu.3G030200.1.p",
             "Urofu.8G010200.1","Sevir.8G008100.1","Chala.12G017400.1",
             "Pahal.3G011900.1.p","Pahal.8G013000.1","Pavag08G014600.1",
             "Sevir.7G316501.1.p","GRMZM2G131516","GRMZM2G015080_P01","Pavag05G020900.1.p")

# C3 species columns
C3_cols <- c("Bradi4g44093.1","LOC_Os11g03110.1","Chala.09G176900.1.p","LOC_Os12g02870.1","OEL13314.1")     


# Dataset 1: Position, Consensus, and C4 species
df_C4 <- df_consensus %>%
  select(Position, Consensus, all_of(C4_cols))

# Dataset 2: Position, Consensus, and C3 species
df_C3 <- df_consensus %>%
  select(Position, Consensus, all_of(C3_cols))


df_C3_clean <- df_C3 %>%
  filter(if_all(all_of(C3_cols), ~ !is.na(.) & . != "-"))

C3_conserved_positions <- df_C3_clean %>%
  # keep only rows where all C3 columns are identical
  filter(apply(select(., all_of(C3_cols)), 1, function(x) length(unique(x)) == 1)) %>%
  pull(Position)

df_C4_conserved_from_C3 <- df_C4 %>%
  filter(Position %in% C3_conserved_positions)

aa_cols <- setdiff(colnames(df_C4_conserved_from_C3), c("Position", "Consensus"))

n_species <- length(aa_cols)

df_consensus_mean <- df_C4_conserved_from_C3 %>%
  rowwise() %>%
  mutate(
    # collect residues for this row
    residues = list(c_across(all_of(aa_cols))),
    
    # count consensus matches
    Consensus_Matches = sum(unlist(residues) == Consensus, na.rm = TRUE),
    
    # count gaps ("-" or "")
    Gaps = sum(unlist(residues) %in% c("-", ""), na.rm = TRUE),
    
    # compute non-gap sequences
    NonGap_Count = n_species - Gaps,
    
    # compute mean consensus score
    Consensus_Mean = ifelse(
      NonGap_Count > 0,
      Consensus_Matches / NonGap_Count,
      NA_real_
    )
  ) %>%
  ungroup()  

df_consensus_mean <- df_consensus_mean %>%
  mutate(
    LowConsensus = Consensus_Mean < 0.4
  )
library(ggplot2)

ggplot(df_consensus_mean, aes(x = Position, y = Consensus_Mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  
  # Normal points
  geom_point(data = subset(df_consensus_mean, !LowConsensus),
             color = "darkblue", size = 2) +
  
  # Highlighted low-consensus positions
  geom_point(data = subset(df_consensus_mean, LowConsensus),
             color = "red", size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Consensus Mean Across Alignment (positions < 0.4 highlighted)",
    x = "Alignment Position",
    y = "Mean Consensus"
  )

ggplot(df_row_scores, aes(x = Position, y = RowSD)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Standard Deviation of Substitution Scores per Alignment Position",
    x = "Alignment Position",
    y = "Row SD (BLOSUM-based)"
  )
aa_cols <- C4_cols
df_row_scores <- df_consensus_mean %>%
  rowwise() %>%
  mutate(
    RowScore = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely from scoring
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Compute substitution scores only for valid AA
      scores <- sapply(valid_residues, function(aa) {
        if (aa == consensus_aa) {
          0
        } else if (aa %in% colnames(blos)) {
          blos[consensus_aa, aa]
        } else {
          NA
        }
      })
      
      # SUM the scores
      sum(scores, na.rm = TRUE)
    },
    
    RowMean = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Number of valid residues (denominator for mean)
      n_valid <- length(valid_residues)
      
      if (n_valid == 0) {
        NA   # no valid AAs at this position
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        sum(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowMean)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "black", size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution BLOSUM62 Score Across Alignment",
    x = "Alignment Position",
    y = "Row BLOSUM62 Substitution Score"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean < 0.4,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_row_scores$Consensus_Mean=1-df_row_scores$Consensus_Mean

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean > 0.75,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )

df_tmp <- df_row_scores %>%
  mutate(
    Low_Frequency = Consensus_Mean > 0.75,
    Low_Score     = RowMean < -0.5
  )
names(df_tmp)
df_tmp2 <- df_tmp
colnames(df_tmp2)[20]=c("Mean substitution frequency")
colnames(df_tmp2)[23]=c("Mean substitution score")
colnames(df_tmp2)[24]=c("Mean substitution standard deviation")

df_long <- df_tmp2 %>%
  select(
    Position,
    `Mean substitution frequency`,
    `Mean substitution score`,
    `Mean substitution standard deviation`,
    Low_Frequency,
    Low_Score
  ) %>%
  pivot_longer(
    cols = c(
      `Mean substitution frequency`,
      `Mean substitution score`,
      `Mean substitution standard deviation`
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = case_when(
      Metric == "Mean substitution frequency"          ~ Low_Frequency,
      Metric == "Mean substitution score"              ~ Low_Score,
      Metric == "Mean substitution standard deviation" ~ FALSE,
      TRUE ~ FALSE
    )
  )

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

msa_seq <- "MGSSSVLLFPSSSSAAPSSHSHATASSHSLLPPLPCDHVLIHYIHQLDEQ-EAATMVRKRPAPDMDLPPPRRHVTGDLSDVTAAAAAGGPGPSSASAQLPALPTQLPAFQHAAEVDVPPAHGGEAPASTTAWVDGIIRDIIGSSGGAVSITQLIHNVREIIHPCNPGLASLLELRLRSLLAADPAPLPQQ-QRALLHGAPAGLALPLPPPLPDKRRHEPAEPHPAPQSPKVPTAEETAA--ASAAAAKERKEVQRRKQRDEEGLHLLTLLLQCAEAVNADNLDDAHQTLLEIAELATPFGTSTQRVAAYFAEAMSARVVSSCLGLYAPLPPGSPAAARLHGRVAAAFQVFNGISPFVKFSHFTANQAIQEAFEREERVHIIDLDIMQGLQWPGLFHILASRPGGPPRVRLTGLGASMEALEATGKRLSDFADTLGLPFEFCAVDEKVGNVDPQKLGVTRREAVAVHWLHHSLYDVTGSDSNTLRLIQRLAPKVVTMVEQDLSQSGSFLARFVDAIHYYSALFDSLDASYGEDSPERHVVEQQLLAREIRNVLAVGGPARAGGARFGSWREELARSGFRAASLAGGAAAQASLLLGMFPSDGYTLVEEKGALRLGWKDLCLLTASAWRPVQTPCR"
maize_seq <- "MGSSSVLLFPSSSSAAPSAPHSFPHSHATAIASSHSLLPPLPCSNPPPPLSSQDHVLIHYIHQLDEQEAATMVRKRPAPDMDLPPPRRHVTGDLSDVTAAAAAGGGPGAPSSASAQLPALPTQLHQLPPAFQHHAAEVDVPPQPHPPAHSQAGGEAPASTTAWVDGIIRDIIGSSGGGAVSITQLIHNVREIIHPCNPGLASLLELRLRSLLAADPAPLPQQQRALLHGAPAAAAGLALPLPPPLPDKRRHEPAPRCQQQQQEEPHPAPQSPKVPTAEETAAASAAAAKERKEVQRRKQRDEEGLHLLTLLLQCAEAVNADNLDDAHQTLLEIAELATPFGTSTQRVAAYFAEAMSARVVSSCLGLYAPLPPGSPAAARLHGRVAAAFQVFNGISPFVKFSHFTANQAIQEAFEREERVHIIDLDIMQGLQWPGLFHILASRPGGPPRVRLTGLGASMEALEATGKRLSDFADTLGLPFEFCAVDEKVGNVDPQKLGVTRREAVAVHWLHHSLYDVTGSDSNTLRLIQRLAPKVVTMVEQDLSQSGSFLARFVDAIHYYSALFDSLDASYGEDSPERHVVEQQLLAREIRNVLAVGGPARAGAGGARFGSWREELARSGFRAASLAGGAAAQASLLLGMFPSDGYTLVEEKGALRLGWKDLCLLTASAWRPVQTPPCR"
mapping_df_2 <- map_msa_to_maize(msa_seq, maize_seq)

domains <- tribble(
  ~Domain, ~Start, ~End,
  "GRAS domain",	298,	630,
  #"bHLH-TF_ACT-like_plant",	165,	241,
  #  "Zinc finger SWIM-type profile",	590,	626,
)
domains

domain_ranges_msa_2 <- domains %>%
  rowwise() %>%
  mutate(
    msa_start = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == Start)[1] ],
    msa_end   = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == End)[1] ]
  ) %>%
  ungroup()


df_score <- df_long %>% 
  filter(Metric == "Mean substitution score")

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

df_long_filtered <- df_long %>%
  filter(Metric != "Mean substitution standard deviation")

df_sd_fixed <- df_sd %>%
  mutate(Metric = "Mean substitution score")

maize_labels <- mapping_df_2$maize_position
names(maize_labels) <- mapping_df_2$msa_position
df_long_filtered <- df_long_filtered %>%
  left_join(mapping_df_2,
            by = c("Position" = "msa_position"))


SCR1 <- ggplot(
  df_long,
  aes(x = Position, y = Value, fill = Metric)
) +
  
  ## --- DOMAIN REGIONS (background) ---
  geom_rect(
    data = domain_ranges_msa_2,
    aes(
      xmin = msa_start - 0.5,
      xmax = msa_end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = Domain
    ),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  
  ## --- MAIN BARS (THICKER) ---
  geom_col(
    width = 1,
    alpha = 0.6,
    position = "identity"
    ) +
  ## maize positions
  scale_x_continuous(
    breaks = mapping_df_2$msa_position,
    labels = maize_labels
  ) +
  ## --- HIGHLIGHTED BARS ---
  geom_point(
    data = df_long_filtered %>% filter(Highlight),
    aes(x = Position, y = Value),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = df_long_filtered %>% filter(Highlight),
    aes(label = maize_position),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  ## --- COLORS ---
  scale_fill_manual(
    values = c(
      "Mean substitution frequency" = "steelblue",
      "Mean substitution score"     = "darkred",
      "GRAS domain"     = "red",
      "bHLH-TF_ACT-like_plant"= "orange" 
    )
  ) +
  
  ## --- AXES ---
  #coord_cartesian(ylim = c(-2, 2), expand = FALSE) +
  
  ## --- LABELS & THEME ---
  theme_minimal(base_size = 14) +
  labs(
    title = "SCR1 & SCR2",
    x = "Alignment Position",
    y = "Value",
    fill = "Metric / Domain"
  ) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
SCR1

#shr1
df=read.csv("MSA_positions_table_GRMZM2G064638.csv", header = T)
# Compute consensus across each row
df_consensus <- df%>%
  rowwise() %>%
  mutate(
    Consensus = {
      # Extract all residues in this row except Position column
      residues <- c_across(-Position)
      
      # Remove gaps or missing characters
      residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # If all sequences have gap, return "-"
      if (length(residues) == 0) {
        "-"
      } else {
        # Find most frequent residue
        tbl <- table(residues)
        names(tbl)[which.max(tbl)]
      }
    }
  ) %>%
  ungroup()

# View result
head(df_consensus)
colnames(df_consensus)

# C4 species columns
C4_cols <- c("ELECO.r07.2AG0108250.1","ELECO.r07.2BG0161890.1","GRMZM2G064638_P01",
             "Pavag04G110200.1","GRMZM2G085751","Urofu.1G118300.1","Pahal.1G118900.1","Sevir.1G013100.1")

# C3 species columns
C3_cols <- c("Bradi3g09670.1","LOC_Os02g15760.1","Chala.11G216500.1","OEL12547.1")     


# Dataset 1: Position, Consensus, and C4 species
df_C4 <- df_consensus %>%
  select(Position, Consensus, all_of(C4_cols))

# Dataset 2: Position, Consensus, and C3 species
df_C3 <- df_consensus %>%
  select(Position, Consensus, all_of(C3_cols))


df_C3_clean <- df_C3 %>%
  filter(if_all(all_of(C3_cols), ~ !is.na(.) & . != "-"))

C3_conserved_positions <- df_C3_clean %>%
  # keep only rows where all C3 columns are identical
  filter(apply(select(., all_of(C3_cols)), 1, function(x) length(unique(x)) == 1)) %>%
  pull(Position)

df_C4_conserved_from_C3 <- df_C4 %>%
  filter(Position %in% C3_conserved_positions)

aa_cols <- setdiff(colnames(df_C4_conserved_from_C3), c("Position", "Consensus"))

n_species <- length(aa_cols)

df_consensus_mean <- df_C4_conserved_from_C3 %>%
  rowwise() %>%
  mutate(
    # collect residues for this row
    residues = list(c_across(all_of(aa_cols))),
    
    # count consensus matches
    Consensus_Matches = sum(unlist(residues) == Consensus, na.rm = TRUE),
    
    # count gaps ("-" or "")
    Gaps = sum(unlist(residues) %in% c("-", ""), na.rm = TRUE),
    
    # compute non-gap sequences
    NonGap_Count = n_species - Gaps,
    
    # compute mean consensus score
    Consensus_Mean = ifelse(
      NonGap_Count > 0,
      Consensus_Matches / NonGap_Count,
      NA_real_
    )
  ) %>%
  ungroup()  

df_consensus_mean <- df_consensus_mean %>%
  mutate(
    LowConsensus = Consensus_Mean < 0.4
  )
library(ggplot2)

ggplot(df_consensus_mean, aes(x = Position, y = Consensus_Mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  
  # Normal points
  geom_point(data = subset(df_consensus_mean, !LowConsensus),
             color = "darkblue", size = 2) +
  
  # Highlighted low-consensus positions
  geom_point(data = subset(df_consensus_mean, LowConsensus),
             color = "red", size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Consensus Mean Across Alignment (positions < 0.4 highlighted)",
    x = "Alignment Position",
    y = "Mean Consensus"
  )

ggplot(df_row_scores, aes(x = Position, y = RowSD)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Standard Deviation of Substitution Scores per Alignment Position",
    x = "Alignment Position",
    y = "Row SD (BLOSUM-based)"
  )
aa_cols <- C4_cols
df_row_scores <- df_consensus_mean %>%
  rowwise() %>%
  mutate(
    RowScore = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely from scoring
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Compute substitution scores only for valid AA
      scores <- sapply(valid_residues, function(aa) {
        if (aa == consensus_aa) {
          0
        } else if (aa %in% colnames(blos)) {
          blos[consensus_aa, aa]
        } else {
          NA
        }
      })
      
      # SUM the scores
      sum(scores, na.rm = TRUE)
    },
    
    RowMean = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Number of valid residues (denominator for mean)
      n_valid <- length(valid_residues)
      
      if (n_valid == 0) {
        NA   # no valid AAs at this position
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        sum(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowMean)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "black", size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution BLOSUM62 Score Across Alignment",
    x = "Alignment Position",
    y = "Row BLOSUM62 Substitution Score"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean < 0.4,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_row_scores$Consensus_Mean=1-df_row_scores$Consensus_Mean

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean > 0.75,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )

df_tmp <- df_row_scores %>%
  mutate(
    Low_Frequency = Consensus_Mean > 0.75,
    Low_Score     = RowMean < -0.5
  )
names(df_tmp)
df_tmp2 <- df_tmp
colnames(df_tmp2)[15]=c("Mean substitution frequency")
colnames(df_tmp2)[18]=c("Mean substitution score")
colnames(df_tmp2)[19]=c("Mean substitution standard deviation")

df_long <- df_tmp2 %>%
  select(
    Position,
    `Mean substitution frequency`,
    `Mean substitution score`,
    `Mean substitution standard deviation`,
    Low_Frequency,
    Low_Score
  ) %>%
  pivot_longer(
    cols = c(
      `Mean substitution frequency`,
      `Mean substitution score`,
      `Mean substitution standard deviation`
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = case_when(
      Metric == "Mean substitution frequency"          ~ Low_Frequency,
      Metric == "Mean substitution score"              ~ Low_Score,
      Metric == "Mean substitution standard deviation" ~ FALSE,
      TRUE ~ FALSE
    )
  )

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

msa_seq <- "MDTLFRLVSLQASEQ-QQQQSASYNSRSTTSSGSRSSSHQTNASYNYYYHSNSSGGGGGQYYYGQQQYLEPYQEECGNTHHLYMDEDFSSSSSSRQHFHSHGAVVQPPTSSTAPTPSLSTSSTAAGAHALFEAADLSFPPDLNLDFSSPASSSGGGAASSAAVGGGGGGRWASQLLLECARAVAGRDSQRVQQLMWMLNELASPYGDVEQKLASYFLQGLFARLTASGPRTLRTLAAASDRNTSFDSTRRTALRFQELSPWSSFGHVAANGAILESFLEAAAASPEPQRLHILDLSNTFCTQWPTLLEALATRSADDTPHLSITTVVSSSAPTAAVQRVMREIGQRMEKFARLMGVPFSFRAVHHAGDLAGLDLDALDLRDGGATTALAINCVNSLRGVVGGARRRDAFAASLRRLDPRVVTVVEEEADLVAFDPGAPESGDTEAAFLKVFGEGLRFFSAYMDSLEESFPKTSNERLALERGAGRAIVDLVSCPASESMERRETAASWARRMRSSGFSPVAFSEDVADDVRSLLRRYREGWSMRDAGLDDSAAGAGVFLAWKEQPLVWASAWRP"
maize_seq <- "MDTLFRLVSLQASEQQQQQSASYNSRSTTSSGSRSSSHQTNASYNYYYHSNSSGGGGGQYYYGQQHPHQHQHQQYYLEPYQQEECGNTHHLYMDEDFSSSSSSRQHFHSHGAVVQPPTSSTATPTAPTPSLSTSSTAAGAAHALFEAADLSFPPDLNLDFSSPASSSGGGAASSAAVGGGGGGRWASQLLLECARAVAGRDSQRVQQLMWMLNELASPYGDVEQKLASYFLQGLFARLTASGPRTLRTLAAASDRNTSFDSTRRTALRFQELSPWSSFGHVAANGAILESFLEAAAASPEPQRLHILDLSNTFCTQWPTLLEALATRSADDTPHLSITTVVSSAPSAPTAAVQRVMREIGQRMEKFARLMGVPFSFRAVHHAGDLAGLDLDALDLRDGGATTALAINCVNSLRGVVPGGARRRDAFAASLRRLDPRVVTVVEEEADLVAFDPGAPEESGDTEAAFLKVFGEGLRFFSAYMDSLEESFPKTSNERLALERGAGRAIVDLVSCPASESMERRETAASWARRMRSSGFSPVAFSEDVADDVRSLLRRYREGWSMRDAGLDDSAAGAGVFLAWKEQPLVWASAWRP"
mapping_df_2 <- map_msa_to_maize(msa_seq, maize_seq)

domains <- tribble(
  ~Domain, ~Start, ~End,
  "GRAS domain",	177,	570,
  #"bHLH-TF_ACT-like_plant",	1,	241,
  #  "Zinc finger SWIM-type profile",	590,	626,
)


domain_ranges_msa_2 <- domains %>%
  rowwise() %>%
  mutate(
    msa_start = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == Start)[1] ],
    msa_end   = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == End)[1] ]
  ) %>%
  ungroup()


df_score <- df_long %>% 
  filter(Metric == "Mean substitution score")

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

df_long_filtered <- df_long %>%
  filter(Metric != "Mean substitution standard deviation")

df_sd_fixed <- df_sd %>%
  mutate(Metric = "Mean substitution score")

maize_labels <- mapping_df_2$maize_position
names(maize_labels) <- mapping_df_2$msa_position
df_long_filtered <- df_long_filtered %>%
  left_join(mapping_df_2,
            by = c("Position" = "msa_position"))


SHR1 <- ggplot(
  df_long,
  aes(x = Position, y = Value, fill = Metric)
) +
  
  ## --- DOMAIN REGIONS (background) ---
  geom_rect(
    data = domain_ranges_msa_2,
    aes(
      xmin = msa_start - 0.5,
      xmax = msa_end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = Domain
    ),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  
  ## --- MAIN BARS (THICKER) ---
  geom_col(
    width = 1,
    alpha = 0.6,
    position = "identity"
    ) +
  ## maize positions
  scale_x_continuous(
    breaks = mapping_df_2$msa_position,
    labels = maize_labels
  ) +
  ## --- HIGHLIGHTED BARS ---
  geom_point(
    data = df_long_filtered %>% filter(Highlight),
    aes(x = Position, y = Value),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = df_long_filtered %>% filter(Highlight),
    aes(label = maize_position),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  ## --- COLORS ---
  scale_fill_manual(
    values = c(
      "Mean substitution frequency" = "steelblue",
      "Mean substitution score"     = "darkred",
      "GRAS domain"     = "red" 
    )
  ) +
  
  ## --- AXES ---
  #coord_cartesian(ylim = c(-2, 2), expand = FALSE) +
  
  ## --- LABELS & THEME ---
  theme_minimal(base_size = 14) +
  labs(
    title = "SHR1 & SHR2",
    x = "Alignment Position",
    y = "Value",
    fill = "Metric / Domain"
  ) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
SHR1

#DOF11
df=read.csv("MSA_positions_table_GRMZM2G123900.csv", header = T)
# Compute consensus across each row
df_consensus <- df%>%
  rowwise() %>%
  mutate(
    Consensus = {
      # Extract all residues in this row except Position column
      residues <- c_across(-Position)
      
      # Remove gaps or missing characters
      residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # If all sequences have gap, return "-"
      if (length(residues) == 0) {
        "-"
      } else {
        # Find most frequent residue
        tbl <- table(residues)
        names(tbl)[which.max(tbl)]
      }
    }
  ) %>%
  ungroup()

# View result
head(df_consensus)
colnames(df_consensus)

# C4 species columns
C4_cols <- c("GRMZM2G123900_P02","GRMZM2G176063","Urofu.3G418800.1",
             "Pavag08G112700.1","Pahal.3G468100.1","Sevir.3G338400.1",
             "ELECO.r07.5BG0442910.1","ELECO.r07.5AG0394670.1")

# C3 species columns
C3_cols <- c("Bradi4g04260.2","LOC_Os12g38200.1","Chala.12G113700.1","OEL22549.1")     


# Dataset 1: Position, Consensus, and C4 species
df_C4 <- df_consensus %>%
  select(Position, Consensus, all_of(C4_cols))

# Dataset 2: Position, Consensus, and C3 species
df_C3 <- df_consensus %>%
  select(Position, Consensus, all_of(C3_cols))


df_C3_clean <- df_C3 %>%
  filter(if_all(all_of(C3_cols), ~ !is.na(.) & . != "-"))

C3_conserved_positions <- df_C3_clean %>%
  # keep only rows where all C3 columns are identical
  filter(apply(select(., all_of(C3_cols)), 1, function(x) length(unique(x)) == 1)) %>%
  pull(Position)

df_C4_conserved_from_C3 <- df_C4 %>%
  filter(Position %in% C3_conserved_positions)

aa_cols <- setdiff(colnames(df_C4_conserved_from_C3), c("Position", "Consensus"))

n_species <- length(aa_cols)

df_consensus_mean <- df_C4_conserved_from_C3 %>%
  rowwise() %>%
  mutate(
    # collect residues for this row
    residues = list(c_across(all_of(aa_cols))),
    
    # count consensus matches
    Consensus_Matches = sum(unlist(residues) == Consensus, na.rm = TRUE),
    
    # count gaps ("-" or "")
    Gaps = sum(unlist(residues) %in% c("-", ""), na.rm = TRUE),
    
    # compute non-gap sequences
    NonGap_Count = n_species - Gaps,
    
    # compute mean consensus score
    Consensus_Mean = ifelse(
      NonGap_Count > 0,
      Consensus_Matches / NonGap_Count,
      NA_real_
    )
  ) %>%
  ungroup()  

df_consensus_mean <- df_consensus_mean %>%
  mutate(
    LowConsensus = Consensus_Mean < 0.4
  )
library(ggplot2)

ggplot(df_consensus_mean, aes(x = Position, y = Consensus_Mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  
  # Normal points
  geom_point(data = subset(df_consensus_mean, !LowConsensus),
             color = "darkblue", size = 2) +
  
  # Highlighted low-consensus positions
  geom_point(data = subset(df_consensus_mean, LowConsensus),
             color = "red", size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Consensus Mean Across Alignment (positions < 0.4 highlighted)",
    x = "Alignment Position",
    y = "Mean Consensus"
  )

ggplot(df_row_scores, aes(x = Position, y = RowSD)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Standard Deviation of Substitution Scores per Alignment Position",
    x = "Alignment Position",
    y = "Row SD (BLOSUM-based)"
  )
aa_cols <- C4_cols
df_row_scores <- df_consensus_mean %>%
  rowwise() %>%
  mutate(
    RowScore = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely from scoring
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Compute substitution scores only for valid AA
      scores <- sapply(valid_residues, function(aa) {
        if (aa == consensus_aa) {
          0
        } else if (aa %in% colnames(blos)) {
          blos[consensus_aa, aa]
        } else {
          NA
        }
      })
      
      # SUM the scores
      sum(scores, na.rm = TRUE)
    },
    
    RowMean = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # Remove gaps/blanks entirely
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      
      # Number of valid residues (denominator for mean)
      n_valid <- length(valid_residues)
      
      if (n_valid == 0) {
        NA   # no valid AAs at this position
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        sum(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

ggplot(df_row_scores, aes(x = Position, y = RowMean)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "black", size = 1.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Substitution BLOSUM62 Score Across Alignment",
    x = "Alignment Position",
    y = "Row BLOSUM62 Substitution Score"
  )
df_row_scores <- df_row_scores %>%
  rowwise() %>%
  mutate(
    RowSD = {
      consensus_aa <- Consensus
      residues     <- c_across(all_of(aa_cols))
      
      # valid amino acids only
      valid_residues <- residues[residues != "-" & residues != "" & !is.na(residues)]
      n_valid <- length(valid_residues)
      
      # SD cannot be computed if fewer than 2 residues
      if (n_valid < 2) {
        NA
      } else {
        scores <- sapply(valid_residues, function(aa) {
          if (aa == consensus_aa) {
            0
          } else if (aa %in% colnames(blos)) {
            blos[consensus_aa, aa]
          } else {
            NA
          }
        })
        
        # Normalized SD
        sd(scores, na.rm = TRUE) / n_valid
      }
    }
  ) %>%
  ungroup()

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean < 0.4,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )


df_row_scores$Consensus_Mean=1-df_row_scores$Consensus_Mean

df_long <- df_row_scores %>%
  mutate(
    Low_Consensus = Consensus_Mean > 0.75,
    Low_RowMean   = RowMean < -0.5
  ) %>%
  select(Position, Consensus_Mean, RowMean, Low_Consensus, Low_RowMean) %>%
  pivot_longer(
    cols = c(Consensus_Mean, RowMean),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight =
      case_when(
        Metric == "Consensus_Mean" & Low_Consensus ~ TRUE,
        Metric == "RowMean"        & Low_RowMean   ~ TRUE,
        TRUE ~ FALSE
      )
  )

df_tmp <- df_row_scores %>%
  mutate(
    Low_Frequency = Consensus_Mean > 0.75,
    Low_Score     = RowMean < -0.5
  )
names(df_tmp)
df_tmp2 <- df_tmp
colnames(df_tmp2)[15]=c("Mean substitution frequency")
colnames(df_tmp2)[18]=c("Mean substitution score")
colnames(df_tmp2)[19]=c("Mean substitution standard deviation")

df_long <- df_tmp2 %>%
  select(
    Position,
    `Mean substitution frequency`,
    `Mean substitution score`,
    `Mean substitution standard deviation`,
    Low_Frequency,
    Low_Score
  ) %>%
  pivot_longer(
    cols = c(
      `Mean substitution frequency`,
      `Mean substitution score`,
      `Mean substitution standard deviation`
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = case_when(
      Metric == "Mean substitution frequency"          ~ Low_Frequency,
      Metric == "Mean substitution score"              ~ Low_Score,
      Metric == "Mean substitution standard deviation" ~ FALSE,
      TRUE ~ FALSE
    )
  )

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

msa_seq <- "M-FLNIPCFEFQPLYIQMQMQQQSPLQCLLGSGGGSDHHHLMPPPSGLAPLPGGPADTAASGPAGGGSSTSAQAAAGAGAQPRPVVSMAERARLARVPLPEPGTLRCPRCDSTNTKFCYFNNYSLSQPRHFCKACRRYWTRGGALRNVPVGGGCRRNTKRSSKKSSRG-GGAGATAATSSSSTTSTSTTATTTTTTSAAMAAAEAIAGMQAQLPHLGLPPAAAAAALEASLEGYHHYLPLQMQPQFLQQAGLHGYHFADDGTGVLADGFPRGVVASGLLAQLAAVKMEEHGSNGGGAIAAHHQSYWPGSTGGGGGWPVEFLSGFSSSSSGNVL"
maize_seq <- "MFLNIPCFEFQPLLIDSLYIQMQMQMQQQSPLQCLLGSGGGSDHHHLMPPPSGLAPLPGGPADTAASGPAGGGSSTSASVQAAAGAGAGAQPRPVVSMAERARLARVPLPEPGTLRCPRCDSTNTKFCYFNNYSLSQPRHFCKACRRYWTRGGALRNVPVGGGCRRNTKRSSKKSSRGGGAGATAATSSSSTTSTSTTATTTTTTSAAMAAAEAIAGMQAQLPHLGLPPAAAAAALEASLEGYHHYLPLQMQPQFLQQAGLHGYHFADDGTGVLAADGFPRGVVASGLLAQLAAVKMEEHGSNGGGAIAAHHEQQSYWPGSTGGGGGWPVEFLSGFSSSSSGNVL*"
mapping_df_2 <- map_msa_to_maize(msa_seq, maize_seq)

domains <- tribble(
  ~Domain, ~Start, ~End,
  "Zinc finger Dof Domain",	115,	169,
  "consensus disorder 1",	34,	94,
  "consensus disorder 2",	159,	201,
)


domain_ranges_msa_2 <- domains %>%
  rowwise() %>%
  mutate(
    msa_start = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == Start)[1] ],
    msa_end   = mapping_df_2$msa_position[ which(mapping_df_2$maize_position == End)[1] ]
  ) %>%
  ungroup()


df_score <- df_long %>% 
  filter(Metric == "Mean substitution score")

df_sd <- df_long %>% 
  filter(Metric == "Mean substitution standard deviation")

df_long_filtered <- df_long %>%
  filter(Metric != "Mean substitution standard deviation")

df_sd_fixed <- df_sd %>%
  mutate(Metric = "Mean substitution score")

maize_labels <- mapping_df_2$maize_position
names(maize_labels) <- mapping_df_2$msa_position
df_long_filtered <- df_long_filtered %>%
  left_join(mapping_df_2,
            by = c("Position" = "msa_position"))


DOF11 <- ggplot(
  df_long,
  aes(x = Position, y = Value, fill = Metric)
) +
  
  ## --- DOMAIN REGIONS (background) ---
  geom_rect(
    data = domain_ranges_msa_2,
    aes(
      xmin = msa_start - 0.5,
      xmax = msa_end + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = Domain
    ),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  
  ## --- MAIN BARS (THICKER) ---
  geom_col(
    width = 1,
    alpha = 0.6,
    position = "identity"
  ) +
  ## maize positions
  scale_x_continuous(
    breaks = mapping_df_2$msa_position,
    labels = maize_labels
  ) +
  ## --- HIGHLIGHTED BARS ---
  geom_point(
    data = df_long_filtered %>% filter(Highlight),
    aes(x = Position, y = Value),
    color = "red",
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = df_long_filtered %>% filter(Highlight),
    aes(label = maize_position),
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  ## --- COLORS ---
  scale_fill_manual(
    values = c(
      "Mean substitution frequency" = "steelblue",
      "Mean substitution score"     = "darkred",
      "Zinc finger Dof Domain" = "red",
      "consensus disorder 1" = "orange",
      "consensus disorder 2" = "lightblue"    )
  ) +
  
  ## --- AXES ---
  #coord_cartesian(ylim = c(-2, 2), expand = FALSE) +
  
  ## --- LABELS & THEME ---
  theme_minimal(base_size = 14) +
  labs(
    title = "DOF11",
    x = "Alignment Position",
    y = "Value",
    fill = "Metric / Domain"
  ) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
DOF11


library(ggpubr)
ggarrange(
  EREB160, DOF11, bhlh116, mads9, bhlh105, SCR1, cadftr3, thx8, c3h28, SHR1,
  ncol = 2,
  nrow = 5,
  heights = rep(2, 6)
)

plots <- list(EREB160, DOF11, bhlh116, mads9, bhlh105, SCR1, cadftr3, thx8, c3h28, SHR1)
plots <- lapply(plots, function(p) 
  p + coord_cartesian(ylim = c(-2, 1))
)

ggarrange(
  plotlist = plots,
  ncol = 2,
  nrow = 6
)


plots <- list(EREB160, DOF11, bhlh116, mads9, bhlh105, SCR1, cadftr3, thx8, c3h28, SHR1)

# remove empty plots (if any)
plots <- Filter(function(p) nrow(p$data) > 0, plots)

# unify y-axis
plots <- lapply(plots, function(p) 
  p + coord_cartesian(ylim = c(-2, 1), expand = FALSE)
)


# arrange
ggarrange(
  plotlist = plots,
  ncol = 2,
  nrow = 5
)
