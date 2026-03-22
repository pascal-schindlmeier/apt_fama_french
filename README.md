# APT Fama-French 3-Faktoren-Modell — Dokumentation

Dieses Skript (`apt_fama_french_lambda_beta.r`) implementiert die **Fama-MacBeth (1973)**-Methodik zur Schätzung der Marktpreise des Risikos (Lambdas) für das Fama-French 3-Faktoren-Modell auf Basis von 25 europäischen Portfolios. Darauf aufbauend wird die erwartete Rendite einer **Einzelaktie** über die APT-Formel berechnet. Das Skript enthält zusätzlich den **GRS-Test** (Gibbons, Ross & Shanken 1989) und das **Cross-sectional R²** als Modellgüte-Diagnostiken.


## Abschnitt 1 — Variablenglossär (chronologisch)

### Benutzereingaben

| Variable | Typ | Beschreibung |
|---|---|---|
| `aktien_name` | `character` | Name der zu analysierenden Aktie (z. B. `"Lufthansa AG"`). Wird ausschließlich für Beschriftungen in der Konsolenausgabe verwendet. |
| `rendite_csv` | `character` | Dateiname der CSV-Datei mit den monatlichen Aktienrenditen. Erwartet zwei Spalten: `Date` (Format YYYYMM) und `Return` (Rendite in Prozent). |

---

### Interne Variablen der Funktion `read_kf_table`

| Variable | Typ | Beschreibung |
|---|---|---|
| `file_path` | `character` | Übergabeparameter: Pfad zur einzulesenden Kenneth-French-Datei. |
| `skip_lines` | `integer` | Übergabeparameter: Anzahl der zu überspringenden Kopfzeilen vor dem eigentlichen Datenteil. |
| `lines` | `character vector` | Alle Rohtextzeilen der Datei, eingelesen via `readLines()`. |
| `start_idx` | `integer` | Erste Zeile des Datenblocks (= `skip_lines + 1`). |
| `empty_lines` | `integer vector` | Indizes aller leeren Zeilen in `lines`, genutzt zur Erkennung des Datenendes. |
| `end_idx` | `integer` | Letzte Zeile des Datenblocks — die Zeile unmittelbar vor der ersten Leerzeile nach `start_idx`. |
| `data_block` | `character vector` | Ausgeschnittener Rohtextvektor, der nur die eigentlichen Datenzeilen enthält. |
| `df` | `data.frame` | Geparster Datensatz aus `data_block`; Rückgabewert der Funktion. Fehlende Werte werden bei `-99.99` und `-999` gesetzt. |

---

### Importierte Rohdaten

| Variable | Typ | Beschreibung |
|---|---|---|
| `ff_factors` | `data.frame` | Fama-French 3-Faktoren für Europa (monatlich). Spalten: `Date` (YYYYMM, numerisch), `Mkt_RF` (Markt-Überschussrendite), `SMB` (Size-Prämie), `HML` (Value-Prämie), `RF` (risikofreier Zinssatz). Alle Werte in Prozent. Quelle: `Europe_3_Factors.csv`. |
| `portfolios` | `data.frame` | Monatliche Renditen der 25 europäischen Portfolios sortiert nach Size × Book-to-Market (5×5-Raster). Spalten: `Date` sowie je eine Spalte pro Portfolio (z. B. `SMALL LoBM`, `BIG HiBM`). Alle Werte in Prozent. Quelle: `Europe_25_Portfolios_ME_BE-ME.csv`. |
| `aktie_returns` | `data.frame` | Monatliche Renditen der Einzelaktie. Spalten: `Date` (numerisch, YYYYMM) und `Return_Aktie` (Rendite in Prozent). Quelle: Datei aus `rendite_csv`. |

---

### Fama-MacBeth Schritt 1 — Zeitreihenregressionen (Portfolio-Betas)

| Variable | Typ | Beschreibung |
|---|---|---|
| `portfolios_long` | `data.frame` | `portfolios` umgewandelt ins Long-Format. Spalten: `Date`, `Portfolio` (Portfolioname als String), `Return` (Rendite in Prozent). Jede Zeile entspricht einer Beobachtung eines Portfolios in einem Monat. |
| `data_portfolios` | `data.frame` | `portfolios_long` innerer Join mit `ff_factors` auf `Date`. Enthält alle Faktoren je Portfolio und Monat sowie die berechnete `Excess_Return`-Spalte (= `Return − RF`). |
| `first_pass` | `data.frame` | Ergebnismatrix der ersten Stufe. Für jedes der 25 Portfolios enthält eine Zeile: `Portfolio` (Name), `alpha_i` (Intercept), `Beta_Mkt`, `Beta_SMB`, `Beta_HML` — jeweils OLS-Schätzer aus der Zeitreihenregression über den gesamten Stichprobenzeitraum. |

---

### Fama-MacBeth Schritt 2 — Querschnittsregressionen (Lambdas)

| Variable | Typ | Beschreibung |
|---|---|---|
| `data_step2` | `data.frame` | Verknüpfung von `data_portfolios` und `first_pass` (Left Join auf `Portfolio`). Enthält für jede Portal-Monats-Beobachtung sowohl die realisierten Renditen als auch die zugehörigen geschätzten Betas. |
| `second_pass` | `data.frame` | Ergebnisse der monatlichen Querschnittsregressionen. Jede Zeile entspricht einem Monat mit den geschätzten Risikopreisen: `Date`, `lambda_0` (Intercept), `lambda_Mkt`, `lambda_SMB`, `lambda_HML`. |

---

### Lambda-Auswertung

| Variable | Typ | Beschreibung |
|---|---|---|
| `lambda_results` | `data.frame` (1 Zeile) | Zusammenfassung der Lambdas über alle Monate. Enthält die zeitlichen Mittelwerte (`lambda_0_avg`, `lambda_Mkt_avg`, `lambda_SMB_avg`, `lambda_HML_avg`) sowie die zugehörigen p-Werte aus einseitigen t-Tests (`p_val_0`, `p_val_Mkt`, `p_val_SMB`, `p_val_HML`). |
| `lambda_table` | `data.frame` | Formatierte Ausgabetabelle der Lambda-Ergebnisse. Spalten: `Faktor` (Bezeichnung), `Lambda_pm` (Mittelwert in % p. M., auf 4 Stellen gerundet), `p_Wert`, `Signifikanz` (Sternchen-Notation nach konventionellen Schwellenwerten). |

---

### Beta-Schätzung Einzelaktie

| Variable | Typ | Beschreibung |
|---|---|---|
| `data_aktie` | `data.frame` | `aktie_returns` innerer Join mit `ff_factors` auf `Date`, ergänzt um `Excess_Return_Aktie` (= `Return_Aktie − RF`). Zeilen mit fehlenden Werten werden entfernt (`drop_na()`). Definiert den Analysezeitraum der Aktie. |
| `aktie_model` | `lm` | OLS-Regressionsmodell der Aktie: `Excess_Return_Aktie ~ Mkt_RF + SMB + HML`. Schätzt die drei Faktorladungen (Betas) und das Alpha der Aktie. |
| `nw_coef` | `coeftest`-Objekt | Koeffizienten von `aktie_model` mit Newey-West HAC-Standardfehlern (Lag = 4). Korrigiert für heteroskedastische und autokorrelierte Residuen. |
| `beta_table` | `data.frame` | Formatierte Ausgabetabelle der Beta-Schätzung. Spalten: `Parameter` (Alpha/Beta_Mkt/Beta_SMB/Beta_HML), `Schaetzer`, `Std_Fehler` (Newey-West), `t_Wert`, `p_Wert`, `Signifikanz`. |
| `adj_r2` | `numeric` | Bereinigtes Bestimmtheitsmaß (Adjusted R²) des `aktie_model`. Gibt an, welcher Anteil der Varianz der Aktienrenditen durch die drei Faktoren erklärt wird. |

---

### APT-Berechnung

| Variable | Typ | Beschreibung |
|---|---|---|
| `beta_mkt` | `numeric` | Markt-Beta der Aktie, extrahiert aus `coef(aktie_model)["Mkt_RF"]`. Misst die Sensitivität der Aktienrendite gegenüber der Markt-Überschussrendite. |
| `beta_smb` | `numeric` | SMB-Beta der Aktie (Size-Faktorladung), extrahiert aus `coef(aktie_model)["SMB"]`. Positiver Wert weist auf Small-Cap-Exposition hin. |
| `beta_hml` | `numeric` | HML-Beta der Aktie (Value-Faktorladung), extrahiert aus `coef(aktie_model)["HML"]`. Positiver Wert weist auf Value-Exposition hin. |
| `lambda_mkt` | `numeric` | Durchschnittlicher Marktpreis des Marktrisikos (`lambda_Mkt_avg`), übernommen aus `lambda_results`. In % p. M. |
| `lambda_smb` | `numeric` | Durchschnittlicher Marktpreis des Size-Risikos (`lambda_SMB_avg`), übernommen aus `lambda_results`. In % p. M. |
| `lambda_hml` | `numeric` | Durchschnittlicher Marktpreis des Value-Risikos (`lambda_HML_avg`), übernommen aus `lambda_results`. In % p. M. |
| `avg_rf` | `numeric` | Durchschnittlicher risikofreier Zinssatz im Analysezeitraum der Aktie (Mittelwert der `RF`-Spalte aus `data_aktie`). In % p. M. |
| `apt_monthly` | `numeric` | APT-geschätzte erwartete Rendite der Aktie in % **pro Monat**: `avg_rf + beta_mkt·lambda_mkt + beta_smb·lambda_smb + beta_hml·lambda_hml`. |
| `apt_annual` | `numeric` | APT-geschätzte erwartete Rendite der Aktie **pro Jahr** (geometrisch aufgezinst): `(1 + apt_monthly/100)^12 − 1`. Wird als Dezimalzahl ausgegeben (×100 für Prozentwert). |
| `apt_breakdown` | `data.frame` | Aufschlüsselung der APT-Ergebnisse nach Komponenten. Spalten: `Komponente` (Bezeichnung), `Beta` (Faktorladung der Aktie), `Lambda_pm` (Risikopreis in % p. M.), `Beitrag_pm` (Beitrag zur erwarteten Rendite in % p. M.). |

---

## Abschnitt 2 — Funktionen und Vorgänge

### `read_kf_table(file_path, skip_lines)`

**Zweck:** Benutzerdefinierter Parser für Datendateien von Kenneth French's Data Library. Diese Dateien enthalten Textköpfe variabler Länge, mehrere Datentabellen, die durch Leerzeilen getrennt sind, und keine standardisierte CSV-Struktur.

**Ablauf:**
1. `readLines()` liest die gesamte Datei als Zeichenvektor (eine Zeile = ein Element).
2. `start_idx` markiert den Beginn des Datenblocks, indem die angegebene Anzahl an Kopfzeilen übersprungen wird.
3. `which(trimws(lines) == "")` findet alle Leerzeilen — das Ende des ersten Datenblocks nach `start_idx` definiert `end_idx`.
4. Der relevante Textblock wird per Indexbereich `lines[start_idx:end_idx]` ausgeschnitten und zu einem einzigen String zusammengefügt.
5. `read_csv()` parst diesen String als CSV; die Sentinel-Werte `-99.99` und `-999` werden als `NA` behandelt.

---

### Dataimport — Abschnitt 2

**`ff_factors`** wird aus `Europe_3_Factors.csv` mit `skip_lines = 6` geladen. Anschließend wird die erste Spalte in `Date` umbenannt, `Mkt-RF` in `Mkt_RF` (Bindestrich zu Unterstrich für R-Kompatibilität), und alle Spalten werden in numerische Typen konvertiert. `Date` wird zunächst als `character` belassen und dann mit `as.numeric()` gewandelt, um Führungsnullen zu vermeiden.

**`portfolios`** wird aus `Europe_25_Portfolios_ME_BE-ME.csv` mit `skip_lines = 20` geladen (längerer Kopfbereich). Die erste Spalte wird in `Date` umbenannt, der Rest bleibt mit Originalnamen erhalten.

**`aktie_returns`** wird direkt mit `read_csv()` geladen (Standard-CSV, kein Custom-Parser nötig). Die Spalten werden einheitlich in `Date` und `Return_Aktie` umbenannt. `Date` wird explizit als `numeric` gespeichert, damit der spätere Join auf `ff_factors` typkompatibel ist.

---

### Fama-MacBeth Schritt 1 — Zeitreihenregressionen

**Ziel:** Für jedes der 25 Portfolios werden portflio-spezifische Faktorladungen (Betas) geschätzt, die messen, wie sensitiv jedes Portfolio gegenüber den drei Fama-French-Faktoren ist.

**Umformung (`portfolios_long`):** `pivot_longer()` bringt die breite Portfoliomatrix (25 Portfliosspalten) in das Long-Format: jede Zeile repräsentiert eine Portfolio-Monats-Beobachtung. Das ist Voraussetzung für Group-by-Operationen.

**Zusammenführung (`data_portfolios`):** `inner_join(ff_factors, by = "Date")` verknüpft jede Portfoliorendite mit den gleichzeitigen Faktorrenditen. Monate ohne Faktorwerte fallen heraus. `Excess_Return = Return − RF` subtrahiert den risikofreien Zins für die eigentliche Regressionsgleichung.

**Regression (`first_pass`):** `group_by(Portfolio) %>% do(tidy(lm(...)))` schätzt für jedes der 25 Portfolios separat die OLS-Regression:

$$\text{ExcessReturn}_{p,t} = \alpha_p + \beta_{p,\text{Mkt}} \cdot \text{Mkt\_RF}_t + \beta_{p,\text{SMB}} \cdot \text{SMB}_t + \beta_{p,\text{HML}} \cdot \text{HML}_t + \varepsilon_{p,t}$$

`broom::tidy()` wandelt das `lm`-Objekt in ein tidy `data.frame` um. `pivot_wider()` bringt die Ergebnisse ins breite Format (eine Portfolio-Zeile, je eine Spalte pro Koeffizient).

---

### Fama-MacBeth Schritt 2 — Querschnittsregressionen

**Ziel:** Monat für Monat wird eine Querschnittsregression über alle 25 Portfolios geschätzt. Die Betas aus Schritt 1 dienen als Regressoren, die realisierten Überschussrenditen als abhängige Variable. Die resulting Koeffizienten sind die monatlichen **Risikopreise (Lambdas)**.

**Verknüpfung (`data_step2`):** `left_join(first_pass, by = "Portfolio")` hängt die zeitkonstanten Betas an jede Zeile des Panel-Datensatzes.

**Monatliche Regression (`second_pass`):** `group_by(Date) %>% do(tidy(lm(...)))` schätzt jeden Monat:

$$\text{ExcessReturn}_{p,t} = \lambda_{0,t} + \lambda_{\text{Mkt},t} \cdot \hat\beta_{p,\text{Mkt}} + \lambda_{\text{SMB},t} \cdot \hat\beta_{p,\text{SMB}} + \lambda_{\text{HML},t} \cdot \hat\beta_{p,\text{HML}} + \varepsilon_{p,t}$$

Pro Monat gibt es 25 Datenpunkte (die 25 Portfolios) und 4 Schätzparameter. Das Ergebnis ist eine Zeitreihe von Lambdas — für jeden der T Monate ein Schätzwert pro Faktor.

---

### Lambda-Auswertung — Abschnitt 5

**`lambda_results`:** `summarize()` berechnet die zeitlichen Mittelwerte der Lambda-Zeitreihen (das eigentliche Fama-MacBeth-Ergebnis) sowie für jeden Lambda-Vektor einen t-Test (`t.test()`), der prüft, ob der Mittelwert signifikant von null verschieden ist. Die Standardabweichung der Lambda-Zeitreihe wird dabei als Schätzer für die Präzision verwendet (Fama-MacBeth-Standardfehler).

**`lambda_table`:** Reine Formatierungsoperation. `mutate(Signifikanz = case_when(...))` fügt Sternchencodes gemäß den konventionellen Schwellenwerten p < 0,01 / 0,05 / 0,10 hinzu.

---

### Beta-Schätzung Einzelaktie — Abschnitt 6

**Datenvorbereitung (`data_aktie`):** Die Aktienrenditen werden per `inner_join` mit den Faktordaten verknüpft (Schnittmenge der Zeiträume). `drop_na()` entfernt Monate mit fehlenden Werten in irgendeiner Spalte.

**OLS-Regression (`aktie_model`):** Zeitreihenregression der Aktienüberschussrenditen auf die drei Fama-French-Faktoren — strukturell identisch zur Portfolio-Regression in Schritt 1, jedoch für eine Einzelaktie:

$$\text{ExcessReturn}_{\text{Aktie},t} = \alpha + \beta_{\text{Mkt}} \cdot \text{Mkt\_RF}_t + \beta_{\text{SMB}} \cdot \text{SMB}_t + \beta_{\text{HML}} \cdot \text{HML}_t + \varepsilon_t$$

**Newey-West-Korrektur (`nw_coef`):** `coeftest(aktie_model, vcov = NeweyWest(aktie_model, lag = 4))` ersetzt die OLS-Standardfehler durch **HAC-Standardfehler** (Heteroskedasticity and Autocorrelation Consistent). Lag = 4 Monate berücksichtigt mögliche Autokorrelation in monatlichen Aktienrenditen. Dies beeinflusst Standardfehler, t-Werte und p-Werte, **nicht** die Punktschätzer der Betas.

**`adj_r2`:** `summary(aktie_model)$adj.r.squared` extrahiert das bereinigte R² aus dem Standard-OLS-Modell (nicht aus `nw_coef`, da HAC-Korrektur nur die Inferenz betrifft).

---

### APT-Berechnung — Abschnitt 7

**Extraktion der Betas:** `coef(aktie_model)` gibt die OLS-Koeffizientenvektoren zurück. Die drei Faktorladungen werden per Namensindizierung (`["Mkt_RF"]` etc.) herausgezogen.

**APT-Formel:** Kombination der aktienspezifischen Betas (aus Abschnitt 6) mit den portfoliobasierten Lambdas (aus Abschnitt 4/5):

$$E(R_i) = \overline{R_F} + \hat\beta_{\text{Mkt}} \cdot \bar\lambda_{\text{Mkt}} + \hat\beta_{\text{SMB}} \cdot \bar\lambda_{\text{SMB}} + \hat\beta_{\text{HML}} \cdot \bar\lambda_{\text{HML}}$$

Dabei stammen β aus der Zeitreihenregression der Einzelaktie und λ aus der Fama-MacBeth-Querschnittsschätzung der 25 Portfolios. `avg_rf` ist der durchschnittliche risikofreie Zinssatz im aktienspezifischen Beobachtungszeitraum.

**Annualisierung (`apt_annual`):** Geometrische Aufzinsung der monatlichen Rendite: `(1 + apt_monthly/100)^12 − 1`. Division durch 100 konvertiert den in Prozent gespeicherten Wert in eine Dezimalzahl vor der Potenzierung.

**`apt_breakdown`:** `data.frame()`-Konstruktion einer Komponentenzerlegung, die den Beitrag jedes Faktors (`Beta × Lambda`) zur erwarteten Gesamtrendite separat ausweist.

---

## Abschnitt 3 — Datenfluss (Mermaid)

flowchart TD
    %% ── KNOTEN-DEFINITIONEN (Inhaltlich unverändert) ────────────────
    %% Eingabedateien
    F1["📄 Europe_3_Factors.csv\n(skip_lines = 6)"]
    F2["📄 Europe_25_Portfolios_ME_BE-ME.csv\n(skip_lines = 20)"]
    F3["📄 lha_monthly_returns.csv\n(Standard-CSV)"]

    %% Kenneth-French-Parser
    P["read_kf_table()\n• readLines → Kopfzeilen überspringen\n• Leerzeile als Tabellenende\n• read_csv auf Datentextblock"]

    %% Zwischendaten & Portfolios
    ff_factors["ff_factors\nDate, Mkt_RF, SMB, HML, RF\n(monatlich, in %)"]
    portfolios["portfolios\nDate + 25 Portfolio-Spalten\n(monatlich, in %)"]
    aktie_returns["aktie_returns\nDate, Return_Aktie\n(monatlich, in %)"]
    
    portfolios_long["portfolios_long\nDate, Portfolio, Return\n(Long-Format)"]
    data_portfolios["data_portfolios\n+ Excess_Return = Return − RF"]
    first_pass["first_pass\nPortfolio | alpha_i | Beta_Mkt | Beta_SMB | Beta_HML\n(25 Zeilen)"]
    data_step2["data_step2\n(Panel: Renditen + Betas)"]
    second_pass["second_pass\nDate | lambda_0 | lambda_Mkt | lambda_SMB | lambda_HML\n(T Zeilen = T Monate)"]
    lambda_results["lambda_results\nlambda_X_avg + p_val_X\n(1 Zeile)"]
    lambda_table["📊 lambda_table\nAusgabe: Risikopreise λ\n(Fama-MacBeth-Ergebnis)"]

    %% Zwischendaten Einzelaktie
    data_aktie["data_aktie\nSchnittmenge Aktie ∩ Faktoren"]
    aktie_model["aktie_model\nOLS-Regressionsobjekt"]
    nw_coef["nw_coef\nKoeffizienten mit robusten SE"]
    adj_r2["adj_r2\nErklärungsgüte"]
    beta_table["📊 beta_table\nAusgabe: Aktien-Betas"]

    %% APT Berechnung
    apt_calc(["APT-Berechnung"])
    apt_monthly["apt_monthly \n Erwartete Rendite p.M. (%)"]
    apt_annual["apt_annual \n Erwartete Rendite p.a. (Dez.)"]
    apt_breakdown["📊 apt_breakdown \n Ausgabe: APT-Ergebnis\nmit Faktorzerlegung"]


    %% ── VERBINDUNGEN (Mit verlängerten Pfeilen gegen Überlappung) ──

    %% Parser
    F1 ---> P
    F2 ---> P
    P ---> ff_factors
    P ---> portfolios
    F3 --->|"read_csv + rename"| aktie_returns

    %% Schritt 1: Long-Format & Join
    portfolios --->|"pivot_longer()"| portfolios_long
    portfolios_long --->|"inner_join(by='Date')"| data_portfolios
    ff_factors --->|"inner_join(by='Date')"| data_portfolios

    %% Schritt 1: Zeitreihenregressionen
    data_portfolios ---->|"group_by(Portfolio)\nOLS: ExcessReturn ~ Mkt_RF + SMB + HML\n(je Portfolio über alle t)"| first_pass

    %% Schritt 2: Querschnittsregressionen
    data_portfolios --->|"left_join(first_pass,\nby='Portfolio')"| data_step2
    first_pass --->|"left_join"| data_step2
    data_step2 ---->|"group_by(Date)\nOLS: ExcessReturn ~ Beta_Mkt + Beta_SMB + Beta_HML\n(je Monat über 25 Portfolios)"| second_pass

    %% Schritt 3: Lambda-Auswertung
    second_pass --->|"summarize(mean + t.test)\nüber alle Monate"| lambda_results
    lambda_results --->|"Formatierung + Signifikanz-Sternchen"| lambda_table

    %% Schritt 4: Beta-Schätzung Einzelaktie
    aktie_returns ---->|"inner_join(by='Date')\n+ Excess_Return_Aktie = Return_Aktie − RF\n+ drop_na()"| data_aktie
    ff_factors --->|"inner_join(by='Date')"| data_aktie

    data_aktie --->|"lm(): ExcessReturn_Aktie ~ Mkt_RF + SMB + HML"| aktie_model
    aktie_model --->|"coeftest + NeweyWest(lag=4)\nHAC-Standardfehler"| nw_coef
    aktie_model --->|"summary()$adj.r.squared"| adj_r2
    nw_coef --->|"Formatierung + Signifikanz-Sternchen"| beta_table

    %% Schritt 5: APT-Berechnung
    aktie_model ---->|"coef(): Beta_Mkt, Beta_SMB, Beta_HML"| apt_calc
    lambda_results ---->|"lambda_Mkt_avg, lambda_SMB_avg, lambda_HML_avg"| apt_calc
    data_aktie ---->|"mean(RF) → avg_rf"| apt_calc

    apt_calc --->|"avg_rf + Σ(β·λ)"| apt_monthly
    apt_calc --->|"(1 + apt_monthly/100)^12 − 1"| apt_annual
    apt_calc --->|"Komponentenzerlegung"| apt_breakdown


    %% ── STYLING ──────────────────────────────────────────────────
    classDef input    fill:#dbeafe,stroke:#3b82f6,color:#1e3a5f
    classDef data     fill:#f0fdf4,stroke:#22c55e,color:#14532d
    classDef model    fill:#fef9c3,stroke:#eab308,color:#713f12
    classDef output   fill:#fce7f3,stroke:#ec4899,color:#831843
    classDef calc     fill:#ede9fe,stroke:#7c3aed,color:#3b0764

    class F1,F2,F3 input
    class ff_factors,portfolios,aktie_returns,portfolios_long,data_portfolios,data_step2,data_aktie data
    class first_pass,second_pass,lambda_results,aktie_model,nw_coef,adj_r2 model
    class lambda_table,beta_table,apt_breakdown output
    class apt_calc,apt_monthly,apt_annual calc