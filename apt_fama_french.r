# ==============================================================================
# EMPIRISCHE KAPITALMARKTANALYSE: Arbitrage Pricing Theory (APT)
# Implementierung des Fama-French 3-Faktoren-Modells via Fama-MacBeth
# ==============================================================================

# 1. PAKETE LADEN
library(tidyverse)
library(broom)
library(lmtest)

# ==============================================================================
# 2. DATEN IMPORTIEREN UND BEREINIGEN (KENNETH FRENCH PARSER)
# Problem: Die CSVs enthalten mehrere Tabellen untereinander. 
# Lösung: Wir lesen nur den Text bis zur ersten leeren Zeile ein.
# ==============================================================================

read_kf_table <- function(file_path, skip_lines) {
  # Lese die Datei zeilenweise als Text ein
  lines <- readLines(file_path, warn = FALSE)
  
  # Start des relevanten Datenblocks (Header)
  start_idx <- skip_lines + 1
  
  # Finde die erste leere Zeile nach dem Start
  empty_lines <- which(trimws(lines) == "")
  end_idx <- min(empty_lines[empty_lines > start_idx]) - 1
  
  # Schneide genau diesen Block aus und wandle ihn in einen Dataframe um
  data_block <- lines[start_idx:end_idx]
  df <- read_csv(paste(data_block, collapse = "\n"), 
                 show_col_types = FALSE, 
                 na = c("-99.99", "-999"))
  return(df)
}

print("Lese saubere Daten ein (nur die Value-Weighted Returns)...")

# --- A) Fama-French 3 Faktoren (Europa) ---
# Header in Zeile 7 -> wir überspringen 6
ff_factors <- read_kf_table("Europe_3_Factors.csv", skip_lines = 6) %>%
  rename(Date = 1, Mkt_RF = `Mkt-RF`) %>%
  mutate(Date = as.character(Date)) %>%
  mutate(across(everything(), as.numeric))

# --- B) 25 Europäische Portfolios (Size & Book-to-Market) ---
# Header in Zeile 21 -> wir überspringen 20
portfolios <- read_kf_table("Europe_25_Portfolios_ME_BE-ME.csv", skip_lines = 20) %>%
  rename(Date = 1) %>%
  mutate(Date = as.character(Date)) %>%
  mutate(across(everything(), as.numeric))


# ==============================================================================
# 3. DATENTRANSFORMATION (Tidy Format & Excess Returns)
# ==============================================================================

# Portfolios vom "breiten" Format ins "lange" Format umwandeln
portfolios_long <- portfolios %>%
  pivot_longer(
    cols = -Date, 
    names_to = "Portfolio", 
    values_to = "Return"
  )

# Zusammenführen von Portfolios und Faktoren
data_merged <- portfolios_long %>%
  inner_join(ff_factors, by = "Date") %>%
  # WICHTIG: APT und CAPM basieren auf ÜBERRENDITEN (Excess Returns).
  mutate(Excess_Return = Return - RF)


# ==============================================================================
# 4. FAMA-MACBETH SCHRITT 1: TIME-SERIES REGRESSION (Zeitreihenregression)
# ==============================================================================

print("Führe Schritt 1 aus: Schätzung der Portfolio-Betas...")

first_pass <- data_merged %>%
  group_by(Portfolio) %>%
  do(tidy(lm(Excess_Return ~ Mkt_RF + SMB + HML, data = .))) %>%
  select(Portfolio, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(
    alpha_i = `(Intercept)`,
    Beta_Mkt = Mkt_RF,
    Beta_SMB = SMB,
    Beta_HML = HML
  )


# ==============================================================================
# 5. FAMA-MACBETH SCHRITT 2: CROSS-SECTIONAL REGRESSION (Querschnittsregression)
# ==============================================================================

print("Führe Schritt 2 aus: Monatliche Querschnittsregressionen...")

data_step2 <- data_merged %>%
  left_join(first_pass, by = "Portfolio")

second_pass <- data_step2 %>%
  group_by(Date) %>%
  do(tidy(lm(Excess_Return ~ Beta_Mkt + Beta_SMB + Beta_HML, data = .))) %>%
  select(Date, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(
    lambda_0 = `(Intercept)`, 
    lambda_Mkt = Beta_Mkt,     
    lambda_SMB = Beta_SMB,     
    lambda_HML = Beta_HML      
  )


# ==============================================================================
# 6. AUSWERTUNG & SIGNIFIKANZTESTS
# ==============================================================================

print("Werte Ergebnisse aus...")

fama_macbeth_results <- second_pass %>%
  ungroup() %>%
  summarize(
    Avg_Lambda_0 = mean(lambda_0, na.rm = TRUE),
    Avg_Lambda_Mkt = mean(lambda_Mkt, na.rm = TRUE),
    Avg_Lambda_SMB = mean(lambda_SMB, na.rm = TRUE),
    Avg_Lambda_HML = mean(lambda_HML, na.rm = TRUE),
    
    p_val_0 = t.test(lambda_0)$p.value,
    p_val_Mkt = t.test(lambda_Mkt)$p.value,
    p_val_SMB = t.test(lambda_SMB)$p.value,
    p_val_HML = t.test(lambda_HML)$p.value
  )

# Ergebnistabelle aufbereiten
results_table <- data.frame(
  Faktor = c("Intercept (Lambda_0)", "Markt (Lambda_Mkt)", "Size (Lambda_SMB)", "Value (Lambda_HML)"),
  Durchschnittliche_Praemie_in_Prozent = c(fama_macbeth_results$Avg_Lambda_0, fama_macbeth_results$Avg_Lambda_Mkt, 
                                           fama_macbeth_results$Avg_Lambda_SMB, fama_macbeth_results$Avg_Lambda_HML),
  p_Wert = c(fama_macbeth_results$p_val_0, fama_macbeth_results$p_val_Mkt, 
             fama_macbeth_results$p_val_SMB, fama_macbeth_results$p_val_HML)
)

results_table <- results_table %>%
  mutate(
    Signifikant = ifelse(p_Wert < 0.05, "Ja", "Nein"),
    Signifikanz_Level = case_when(
      p_Wert < 0.01 ~ "***",
      p_Wert < 0.05 ~ "**",
      p_Wert < 0.10 ~ "*",
      TRUE ~ ""
    )
  )

print("=== FINALE ERGEBNISSE FÜR DIE PRÄSENTATION ===")
print(results_table)