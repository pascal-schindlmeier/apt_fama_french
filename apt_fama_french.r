# ==============================================================================
# EMPIRISCHE KAPITALMARKTANALYSE: Arbitrage Pricing Theory (APT)
# Fama-French 3-Faktoren-Modell via Fama-MacBeth + Beta-Schätzung Einzelaktie
# ==============================================================================
#
# METHODIK-ÜBERBLICK:
#   Stufe 1 (Zeitreihe):    Für jedes der 25 Portfolios wird eine OLS-Regression
#                           auf die 3 Fama-French-Faktoren geschätzt → Betas.
#   Stufe 2 (Querschnitt): Jeden Monat: OLS über die 25 Portfolios mit Betas als
#                           Regressoren → monatliche Risikopreise (Lambdas).
#   Stufe 3 (Mittelung):   Zeitlicher Mittelwert der Lambdas = Fama-MacBeth-
#                           Schätzer. t-Test auf Signifikanz (mit na.omit() für
#                           korrekte NA-Behandlung).
#   Diagnostik:            GRS-Test (H0: alle Alpha_i = 0) prüft Modellgüte.
#                           Cross-sectional R² misst Erklärungskraft der Betas.
#   APT-Formel:            E(R_i) = RF + λ_0 + β_Mkt·λ_Mkt + β_SMB·λ_SMB
#                                       + β_HML·λ_HML
#                           λ_0 ist der Intercept der Querschnittsregression.
#                           Er wird explizit aufgenommen, da er empirisch oft
#                           signifikant von null verschieden ist (Misspezifikation
#                           oder Risikoprämien nicht im Modell erfasst).
#
# BEKANNTE EINSCHRÄNKUNGEN (dokumentiert, nicht implementiert):
#   - Full-sample Betas: Betas werden aus dem gesamten Zeitraum geschätzt und
#     in allen Monaten der zweiten Stufe verwendet. Dies erzeugt einen
#     Look-ahead Bias. In der Forschung (inkl. Fama & MacBeth 1973 selbst)
#     ist dies für Risikopreisschätzungen (nicht Out-of-Sample) akzeptiert.
#     Alternative: rollierende Beta-Schätzung (z.B. 60-Monats-Fenster).
#   - Shanken-Korrektur (1992): Die Fama-MacBeth-Standardfehler behandeln die
#     Betas als wahre Werte. Da sie geschätzt sind (EIV-Problem), sind die
#     FM-Standardfehler für Lambdas zu eng. Shanken (1992) liefert den
#     Korrekturfaktor: SE_korr = SE_FM * sqrt(1 + λ' Σ_f⁻¹ λ).
#     Dieser ist > 1 und macht Signifikanzaussagen konservativer. Für
#     explorative Analysen ist der einfache FM-t-Test üblich.
# ==============================================================================

# ==============================================================================
# INPUT: Hier anpassen
# ==============================================================================

aktien_name <- "Lufthansa AG"
rendite_csv <- "lha_monthly_returns.csv"   # Spalten: Date (YYYYMM), Return (%)

# ==============================================================================
# 1. PAKETE LADEN
# ==============================================================================

library(tidyverse)
library(broom)
library(lmtest)
library(sandwich)

# ==============================================================================
# 2. DATEN IMPORTIEREN (KENNETH FRENCH PARSER)
# ==============================================================================

# Kenneth-French-Dateien haben einen Textkopf variabler Länge und mehrere
# Datentabellen, getrennt durch Leerzeilen. Diese Funktion extrahiert den
# ersten Datenblock ab skip_lines und parst ihn als CSV.
# Sentinel-Werte -99.99 und -999 werden als NA behandelt.
read_kf_table <- function(file_path, skip_lines) {
  lines      <- readLines(file_path, warn = FALSE)
  start_idx  <- skip_lines + 1
  empty_lines <- which(trimws(lines) == "")
  # Erstes Leerzeichen NACH dem Datenbeginn markiert das Tabellenende
  end_idx    <- min(empty_lines[empty_lines > start_idx]) - 1
  data_block <- lines[start_idx:end_idx]
  df <- read_csv(paste(data_block, collapse = "\n"),
                 show_col_types = FALSE,
                 na = c("-99.99", "-999"))
  return(df)
}

print("Lese Daten ein...")

# --- A) Fama-French 3 Faktoren (Europa) ---
# Bindestrich in "Mkt-RF" → Unterstrich für R-Kompatibilität.
# Date wird als character gelesen, dann numeric, um Führungsnullen zu vermeiden.
ff_factors <- read_kf_table("Europe_3_Factors.csv", skip_lines = 6) %>%
  rename(Date = 1, Mkt_RF = `Mkt-RF`) %>%
  mutate(Date = as.character(Date)) %>%
  mutate(across(everything(), as.numeric))

# --- B) 25 Europäische Portfolios (Size × Book-to-Market, 5×5-Raster) ---
# skip_lines = 20 überspringt den längeren Textkopf dieser Datei.
portfolios <- read_kf_table("Europe_25_Portfolios_ME_BE-ME.csv", skip_lines = 20) %>%
  rename(Date = 1) %>%
  mutate(Date = as.character(Date)) %>%
  mutate(across(everything(), as.numeric))

# --- C) Einzelaktie (Standard-CSV, kein Custom-Parser nötig) ---
# Date als numeric, damit der spätere Join auf ff_factors typkompatibel ist.
aktie_returns <- read_csv(rendite_csv, show_col_types = FALSE) %>%
  rename(Date = 1, Return_Aktie = 2) %>%
  mutate(Date = as.numeric(Date))


# ==============================================================================
# 3. FAMA-MACBETH SCHRITT 1: ZEITREIHENREGRESSION (25 Portfolios → Betas)
#
# Für jedes der 25 Portfolios wird OLS über die gesamte Zeitreihe geschätzt:
#   Excess_Return_p,t = α_p + β_Mkt·Mkt_RF_t + β_SMB·SMB_t + β_HML·HML_t + ε
#
# Ergebnis: Portfolio-spezifische Betas (β_Mkt, β_SMB, β_HML) und Alphas.
#
# HINWEIS Look-ahead Bias: Die Betas nutzen den gesamten Stichprobenzeitraum,
# werden aber in alle monatlichen Querschnittsregressionen (Schritt 4) eingesetzt
# — auch für frühe Monate, deren Renditen die Beta-Schätzung mitgeformt haben.
# Für die Schätzung von Risikopreisen (nicht Prognose) ist dies in der Literatur
# akzeptiert (vgl. Fama & MacBeth 1973, Cochrane 2005, Kap. 12).
# ==============================================================================

print("Fama-MacBeth Schritt 1: Schätzung der Portfolio-Betas (full-sample OLS)...")

# Breites Format → Long-Format: eine Zeile pro Portfolio × Monat
portfolios_long <- portfolios %>%
  pivot_longer(cols = -Date, names_to = "Portfolio", values_to = "Return")

# Verknüpfung mit Faktordaten; innerer Join → nur Monate mit Faktordaten bleiben
data_portfolios <- portfolios_long %>%
  inner_join(ff_factors, by = "Date") %>%
  mutate(Excess_Return = Return - RF)   # Überschussrendite: r_p - r_f

# group_modify() als moderner Ersatz für das veraltete do() (deprecated seit dplyr 1.0)
first_pass <- data_portfolios %>%
  group_by(Portfolio) %>%
  group_modify(~ tidy(lm(Excess_Return ~ Mkt_RF + SMB + HML, data = .x))) %>%
  select(Portfolio, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(
    alpha_i  = `(Intercept)`,   # Jensen's Alpha: Überrendite gegenüber Modellvorhersage
    Beta_Mkt = Mkt_RF,          # Marktrisiko-Exposition
    Beta_SMB = SMB,             # Size-Exposition (positiv = Small-Cap)
    Beta_HML = HML              # Value-Exposition (positiv = Value-Aktie)
  )


# ==============================================================================
# 4. FAMA-MACBETH SCHRITT 2: QUERSCHNITTSREGRESSION (Betas → Lambdas)
#
# Jeden Monat t wird über die 25 Portfolios geschätzt:
#   Excess_Return_p,t = λ_0,t + λ_Mkt,t·β_Mkt,p + λ_SMB,t·β_SMB,p
#                               + λ_HML,t·β_HML,p + ε_p,t
#
# Die β kommen aus Schritt 1 (zeitkonstant). Pro Monat gibt es 25 Datenpunkte
# und 4 Parameter → 21 Freiheitsgrade. Das Ergebnis ist eine T-lange Zeitreihe
# monatlicher Risikopreisschätzer λ_t für jeden der drei Faktoren plus Intercept.
# ==============================================================================

print("Fama-MacBeth Schritt 2: Monatliche Querschnittsregressionen (25 Portfolios)...")

# Zeitkonstante Betas an jede Zeile des Panel-Datensatzes anhängen
data_step2 <- data_portfolios %>%
  left_join(first_pass, by = "Portfolio")

# group_modify(): je eine OLS-Regression pro Monat über die 25 Portfolios
second_pass <- data_step2 %>%
  group_by(Date) %>%
  group_modify(~ tidy(lm(Excess_Return ~ Beta_Mkt + Beta_SMB + Beta_HML, data = .x))) %>%
  select(Date, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(
    lambda_0   = `(Intercept)`,  # Intercept der Querschnittsregression (≠ 0 → Misspezifikation)
    lambda_Mkt = Beta_Mkt,       # Monatlicher Risikopreis des Marktfaktors
    lambda_SMB = Beta_SMB,       # Monatlicher Risikopreis des Size-Faktors
    lambda_HML = Beta_HML        # Monatlicher Risikopreis des Value-Faktors
  )


# ==============================================================================
# 5. LAMBDA-AUSWERTUNG (Risikopreise & Signifikanz)
#
# Der Fama-MacBeth-Schätzer ist der zeitliche Mittelwert der monatlichen Lambdas:
#   λ̄ = (1/T) Σ λ_t
#
# Signifikanz: t-Test H0: λ̄ = 0 mit FM-Standardfehler = SD(λ_t) / √T.
# na.omit() korrigiert den Bug, dass t.test() bei NA-Werten silent NaN liefert.
#
# Hinweis Shanken-Korrektur (1992): Da die Betas geschätzt (nicht wahr) sind,
# unterschätzt der FM-Standardfehler die wahre Schätzunsicherheit der Lambdas
# (Errors-in-Variables). Shanken (1992) zeigt: SE_korr = SE_FM · √(1 + λ'Σ_f⁻¹λ).
# Der Korrekturfaktor ist > 1 → konservativere Signifikanzaussagen. Für
# explorative Analysen sind die FM-t-Tests dennoch die Standardmethode.
# ==============================================================================

print("Werte Lambda-Schätzungen aus (Fama-MacBeth-Mittelwerte + t-Tests)...")

lambda_results <- second_pass %>%
  ungroup() %>%
  summarize(
    lambda_0_avg   = mean(lambda_0,   na.rm = TRUE),
    lambda_Mkt_avg = mean(lambda_Mkt, na.rm = TRUE),
    lambda_SMB_avg = mean(lambda_SMB, na.rm = TRUE),
    lambda_HML_avg = mean(lambda_HML, na.rm = TRUE),
    # na.omit() ist kritisch: t.test() gibt silent NaN zurück wenn NA-Werte
    # in der Zeitreihe vorhanden sind (z.B. durch degenerierte Monatssregressionen)
    p_val_0        = t.test(na.omit(lambda_0))$p.value,
    p_val_Mkt      = t.test(na.omit(lambda_Mkt))$p.value,
    p_val_SMB      = t.test(na.omit(lambda_SMB))$p.value,
    p_val_HML      = t.test(na.omit(lambda_HML))$p.value
  )

lambda_table <- data.frame(
  Faktor    = c("Lambda_0 (Intercept)", "Lambda_Mkt (Markt)",
                "Lambda_SMB (Size)", "Lambda_HML (Value)"),
  Lambda_pm = round(c(lambda_results$lambda_0_avg, lambda_results$lambda_Mkt_avg,
                      lambda_results$lambda_SMB_avg, lambda_results$lambda_HML_avg), 4),
  p_Wert    = round(c(lambda_results$p_val_0, lambda_results$p_val_Mkt,
                      lambda_results$p_val_SMB, lambda_results$p_val_HML), 4)
) %>%
  mutate(
    Signifikanz = case_when(
      p_Wert < 0.01 ~ "***",
      p_Wert < 0.05 ~ "**",
      p_Wert < 0.10 ~ "*",
      TRUE          ~ ""
    )
  )

print("=== FAMA-MACBETH ERGEBNISSE: RISIKOPREISE (LAMBDAS) ===")
print(lambda_table)
cat("Signifikanzniveaus: *** p<0.01  ** p<0.05  * p<0.10\n")
cat("Standardfehler: Fama-MacBeth (SD(λ_t)/√T) — ohne Shanken-Korrektur\n\n")


# ==============================================================================
# 5a. MODELLGÜTE-DIAGNOSTIKEN
#
# (i)  GRS-Test (Gibbons, Ross & Shanken 1989):
#      Prüft gemeinsam H0: α_1 = α_2 = ... = α_25 = 0 (alle Portfolio-Alphas null).
#      Ein signifikanter GRS-F-Wert zeigt, dass die drei Faktoren die
#      Querschnittsvariation der 25 Portfolio-Renditen NICHT vollständig erklären
#      → das Modell ist misspezifiziert.
#
# (ii) Cross-sectional R²:
#      Regression der zeitlich gemittelten Portfolio-Überschussrenditen auf die
#      mittleren Beta-Kombinationen. Misst, wie gut die geschätzten Betas die
#      durchschnittlichen Renditeunterschiede zwischen den 25 Portfolios erklären.
# ==============================================================================

print("Modellgüte-Diagnostiken: GRS-Test und Cross-sectional R²...")

# --- (i) GRS-Test ---
# Benötigte Größen: Alpha-Vektor (25×1), Faktor-Kovarianzmatrix (3×3),
# Residualkovarianzmatrix (25×25), mittlere Faktorrenditen (3×1).

# Zeitliche Mittelwerte der drei Faktoren (für GRS-Formel)
T_obs    <- length(unique(data_portfolios$Date))
mu_f     <- colMeans(ff_factors[, c("Mkt_RF", "SMB", "HML")], na.rm = TRUE)
Sigma_f  <- cov(ff_factors[, c("Mkt_RF", "SMB", "HML")], use = "complete.obs")

# Alpha-Vektor aus first_pass (25 Portfolio-Alphas)
alphas   <- first_pass$alpha_i

# Residualkovarianzmatrix: OLS-Residuen der 25 Zeitreihenregressionen
resid_matrix <- data_portfolios %>%
  group_by(Portfolio) %>%
  group_modify(~ {
    fit <- lm(Excess_Return ~ Mkt_RF + SMB + HML, data = .x)
    tibble(resid = residuals(fit), Date = .x$Date)
  }) %>%
  ungroup() %>%
  pivot_wider(names_from = Portfolio, values_from = resid) %>%
  select(-Date)

Sigma_eps <- cov(resid_matrix, use = "complete.obs")

# GRS-F-Statistik (Gibbons, Ross & Shanken 1989, Econometrica):
#   F = (T-N-K)/N · (1 + μ_f' Σ_f⁻¹ μ_f)⁻¹ · α' Σ_ε⁻¹ α
#   mit T = Monate, N = 25 Portfolios, K = 3 Faktoren
N_port  <- length(alphas)
K_fact  <- 3
sh2_f   <- as.numeric(t(mu_f) %*% solve(Sigma_f) %*% mu_f)  # Quadratisches Sharpe-Ratio Faktoren

grs_F   <- ((T_obs - N_port - K_fact) / N_port) *
           (1 / (1 + sh2_f)) *
           as.numeric(t(alphas) %*% solve(Sigma_eps) %*% alphas)

grs_p   <- pf(grs_F, df1 = N_port, df2 = T_obs - N_port - K_fact, lower.tail = FALSE)

cat("=== GRS-TEST: H0: Alle 25 Portfolio-Alphas = 0 ===\n")
cat(sprintf("GRS F-Statistik : %.4f\n", grs_F))
cat(sprintf("GRS p-Wert      : %.4f%s\n", grs_p,
            ifelse(grs_p < 0.01, " ***", ifelse(grs_p < 0.05, " **",
                   ifelse(grs_p < 0.10, " *", "")))))
cat(sprintf("Freiheitsgrade  : df1 = %d, df2 = %d\n", N_port, T_obs - N_port - K_fact))
cat("Interpretation: Signifikant → Faktoren erklären Querschnitt NICHT vollständig\n\n")

# --- (ii) Cross-sectional R² ---
# Mittlere Überschussrenditen der 25 Portfolios vs. mittlere Beta-Kombination
mean_returns <- data_portfolios %>%
  group_by(Portfolio) %>%
  summarize(Mean_Excess_Return = mean(Excess_Return, na.rm = TRUE), .groups = "drop")

cs_data <- first_pass %>%
  left_join(mean_returns, by = "Portfolio")

# OLS: Mean_Excess_Return ~ Beta_Mkt + Beta_SMB + Beta_HML (kein Intercept-Zwang)
cs_model  <- lm(Mean_Excess_Return ~ Beta_Mkt + Beta_SMB + Beta_HML, data = cs_data)
cs_r2     <- summary(cs_model)$r.squared
cs_adj_r2 <- summary(cs_model)$adj.r.squared

cat("=== CROSS-SECTIONAL R²: Betas erklären ⌀-Renditeunterschiede der 25 Portfolios ===\n")
cat(sprintf("R²       : %.4f (%.1f%%)\n", cs_r2, cs_r2 * 100))
cat(sprintf("Adj. R²  : %.4f (%.1f%%)\n", cs_adj_r2, cs_adj_r2 * 100))
cat("Interpretation: Anteil der ⌀-Renditevariation zwischen den 25 Portfolios,\n")
cat("                der durch die drei Faktor-Betas erklärt wird.\n\n")


# ==============================================================================
# 6. BETA-SCHÄTZUNG EINZELAKTIE (Zeitreihenregression)
#
# Strukturell identisch zu Schritt 3, jedoch für eine Einzelaktie statt eines
# Portfolios. Faktorladungen (Betas) und Alpha der Aktie werden geschätzt:
#   Excess_Return_Aktie,t = α + β_Mkt·Mkt_RF_t + β_SMB·SMB_t + β_HML·HML_t + ε_t
#
# Standardfehler: Newey-West HAC (lag=4) korrigiert für Heteroskedastizität und
# Autokorrelation in monatlichen Einzelaktienrenditen. Dies betrifft nur Inferenz
# (SE, t, p), NICHT die Punktschätzer der Betas.
# ==============================================================================

cat(sprintf("Schätze Fama-French Betas für %s...\n", aktien_name))

# inner_join: Schnittmenge von Aktien- und Faktorzeitraum
data_aktie <- aktie_returns %>%
  inner_join(ff_factors, by = "Date") %>%
  mutate(Excess_Return_Aktie = Return_Aktie - RF) %>%
  drop_na()   # Entfernt Monate mit fehlenden Werten in einer beliebigen Spalte

cat(sprintf("Analysezeitraum Aktie: %s bis %s (%d Monate)\n",
            min(data_aktie$Date), max(data_aktie$Date), nrow(data_aktie)))

aktie_model <- lm(Excess_Return_Aktie ~ Mkt_RF + SMB + HML, data = data_aktie)

# Newey-West HAC-Standardfehler (Heteroskedasticity and Autocorrelation Consistent)
# lag=4: berücksichtigt Autokorrelation bis zu 4 Monate (übliche Wahl für Monatsdaten)
nw_coef <- coeftest(aktie_model, vcov = NeweyWest(aktie_model, lag = 4))

beta_table <- data.frame(
  Parameter  = c("Alpha", "Beta_Mkt", "Beta_SMB", "Beta_HML"),
  Schaetzer  = round(nw_coef[, "Estimate"], 4),
  Std_Fehler = round(nw_coef[, "Std. Error"], 4),    # Newey-West SE
  t_Wert     = round(nw_coef[, "t value"], 3),
  p_Wert     = round(nw_coef[, "Pr(>|t|)"], 4)
) %>%
  mutate(
    Signifikanz = case_when(
      p_Wert < 0.01 ~ "***",
      p_Wert < 0.05 ~ "**",
      p_Wert < 0.10 ~ "*",
      TRUE          ~ ""
    )
  )

rownames(beta_table) <- NULL
cat(sprintf("\n=== BETA-SCHÄTZUNG: %s ===\n", toupper(aktien_name)))
print(beta_table)

# adj.r.squared aus dem OLS-Basismodell (nicht aus nw_coef, da HAC nur Inferenz betrifft)
adj_r2 <- summary(aktie_model)$adj.r.squared
cat(sprintf("\nAdj. R²: %.4f (%.1f%% der Varianz erklärt)\n", adj_r2, adj_r2 * 100))
cat("Standardfehler: Newey-West HAC (lag=4) — korrigiert für Heteroskedastizität und Autokorrelation\n\n")


# ==============================================================================
# 7. APT-FORMEL
#
# Die APT-Formel setzt die aktienspezifischen Betas (Schritt 6) mit den
# portfoliobasierten Risikopreisen (Schritt 4/5) zusammen:
#
#   E(R_i) = RF + λ_0 + β_Mkt·λ_Mkt + β_SMB·λ_SMB + β_HML·λ_HML
#
# WARUM λ_0 eingeschlossen?
#   Die Querschnittsregression läuft auf Überschussrenditen. Ein von null
#   verschiedener λ_0 signalisiert, dass das Modell nicht alle Risikoprämien
#   erfasst. Wenn λ_0 in lambda_table signifikant ist, MUSS er in die APT-Rendite
#   einbezogen werden — ihn wegzulassen würde die erwartete Rendite systematisch
#   unterschätzen (bei λ_0 > 0) oder überschätzen (bei λ_0 < 0).
#
# β: OLS-Schätzer aus Zeitreihenregression der Einzelaktie (Schritt 6)
# λ: Fama-MacBeth-Mittelwerte aus Querschnittsregression (Schritt 4/5)
# RF: durchschnittlicher risikofreier Zins im Aktien-Analysezeitraum
# ==============================================================================

cat(sprintf("=== APT-SCHÄTZUNG: ERWARTETE RENDITE %s ===\n", toupper(aktien_name)))

# OLS-Betas der Aktie (Punktschätzer, nicht Newey-West-korrigiert)
beta_mkt <- coef(aktie_model)["Mkt_RF"]
beta_smb <- coef(aktie_model)["SMB"]
beta_hml <- coef(aktie_model)["HML"]

# Fama-MacBeth-Risikopreise aus Schritt 5
lambda_0   <- lambda_results$lambda_0_avg    # Intercept-Prämie (Misspezifikation)
lambda_mkt <- lambda_results$lambda_Mkt_avg
lambda_smb <- lambda_results$lambda_SMB_avg
lambda_hml <- lambda_results$lambda_HML_avg

# Durchschnittlicher risikofreier Zins im aktienspezifischen Analysezeitraum
avg_rf <- mean(data_aktie$RF, na.rm = TRUE)

# APT: E(R_i) = RF + λ_0 + β_Mkt·λ_Mkt + β_SMB·λ_SMB + β_HML·λ_HML
# λ_0 wird explizit addiert (war vorher fehlerhaft weggelassen)
apt_monthly <- avg_rf + lambda_0 +
               beta_mkt * lambda_mkt +
               beta_smb * lambda_smb +
               beta_hml * lambda_hml

# Geometrische Aufzinsung: apt_monthly ist in %, daher /100 vor Potenzierung
apt_annual <- (1 + apt_monthly / 100)^12 - 1

apt_breakdown <- data.frame(
  Komponente = c("Risikoloser Zins (RF)",
                 "Intercept-Prämie (λ_0)",
                 "Markt:  β_Mkt · λ_Mkt",
                 "Size:   β_SMB · λ_SMB",
                 "Value:  β_HML · λ_HML"),
  Beta       = c(NA, NA,
                 round(beta_mkt, 4), round(beta_smb, 4), round(beta_hml, 4)),
  Lambda_pm  = c(NA, round(lambda_0, 4),
                 round(lambda_mkt, 4), round(lambda_smb, 4), round(lambda_hml, 4)),
  Beitrag_pm = round(c(avg_rf, lambda_0,
                        beta_mkt * lambda_mkt,
                        beta_smb * lambda_smb,
                        beta_hml * lambda_hml), 4)
)

print(apt_breakdown)

cat(sprintf("\nErwartete Rendite %s (p.M.): %.4f %%\n", aktien_name, apt_monthly))
cat(sprintf("Erwartete Rendite %s (p.a.): %.2f %%\n", aktien_name, apt_annual * 100))
cat(sprintf("\nFormel: E(R_%s) = RF + λ_0 + β_Mkt·λ_Mkt + β_SMB·λ_SMB + β_HML·λ_HML\n", aktien_name))
cat(sprintf("        (λ aus Fama-MacBeth / 25 Portfolios | β aus Zeitreihenregression %s)\n", aktien_name))
cat(sprintf("        λ_0 = %.4f %% p.M. (Intercept-Prämie aus Querschnittsregression)\n",
            lambda_0))
